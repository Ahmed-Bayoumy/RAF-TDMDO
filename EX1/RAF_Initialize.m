%-------------------------------------------------------------------------------------%
%  RAF_MADS                                                                           %
%                                                                                     %
%  A simple matlab version of the Relative Adequacy Framework (RAF) algorithm         %
%  implemented in conjunction with the Mesh Adaptive Direct Search algorithm          %
%  for constrained derivative free optimization.                                      %
%  Version 1.0.0                                                                      %
%                                                                                     %
%  Copyright (C) 2015-2020  Ahmed Bayoumy - McGill University, Montreal               %
%                                                                                     %
%  Author: Ahmed Bayoumy                                                              %
%  email: ahmed.bayoumy@mcgill.ca                                                     %
%                                                                                     %
%  This program is free software: you can redistribute it and/or modify it under the  %
%  terms of the GNU Lesser General Public License as published by the Free Software   %
%  Foundation, either version 3 of the License, or (at your option) any later         %
%  version.                                                                           %
%                                                                                     %
%  This program is distributed in the hope that it will be useful, but WITHOUT ANY    %
%  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A    %
%  PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.   %
%                                                                                     %
%  You should have received a copy of the GNU Lesser General Public License along     %
%  with this program. If not, see <http://www.gnu.org/licenses/>.                     %
%                                                                                     %
%  You can find information on RAF_MADS at                                            %
%  https://github.com/Ahmed-Bayoumy/RAF-TDMDO.git                                     %
%-------------------------------------------------------------------------------------%

%-------------------------------------------------------------------------------------%
%                 /\\\\\\\          /\       /\\\\\\\\                                % 
%                 /\\    /\\       /\ \\     /\\                                      %
%                 /\\    /\\      /\  /\\    /\\                                      %
%                 /\ /\\         /\\   /\\   /\\\\\\                                  %
%                 /\\  /\\      /\\\\\\ /\\  /\\                                      %
%                 /\\    /\\   /\\       /\\ /\\                                      %
%                 /\\      /\\/\\         /\\/\\                                      %
%-------------------------------------------------------------------------------------%
                                    


function [outputArg1,outputArg2] = RAF_Initialize



%% Models definitions

% MF1=[pwd '\MS\MF1.wbpj'];
% MF2=[pwd '\MS\MF2.wbpj'];
% MF3=[pwd '\MS\MF3.wbpj'];
% MF4=[pwd '\MS\MF4.wbpj'];
% 
% MS1=[pwd '\MS\MS1.wbpj'];
% MS2=[pwd '\MS\MS2.wbpj'];
% MS3=[pwd '\MS\MS3.wbpj'];
% MS4=[pwd '\MS\MS4.wbpj'];

MSi=[1 1 1;
    2 1 2;
    3 2 2;
    4 2 3];
MFi=[1 1 2;
    2 1 3;
    3 2 2;
    4 2 1];

MM.MS.APP={'transient', 'static'};
MM.MS.size=[0.01;0.03;0.05];
MM.MF.size=[0.03;0.05];
MM.MF.DT=[0.1;0.01;0.05];

% Predefined cost ratios
MM.MS.CR=[1;0.6;0.2;0.12];
MM.MF.CR=[1;0.73;0.17;0.09];

%% Initializations
disp('========== Initialization ==============');
xlb = [0;0]';
xub = [1;1]';

N=2;
px0 = 6;
LH = lhsdesign(px0,N);
% x0 = ones(px0,1)*xlb0 + LH.*( ones(px0,1)*(xub0-xlb0) );
x0=[0.5 0.5];
OutM=[1 1];
Ef=zeros;
Es=zeros;

Inputfile='MI.wbpjn'; % Standard ANSYS script with tagged design and modeling parameers

for i=1:1:4
    for j=1:1:length(LH)
        Minx = [i 1];
        OutputFile=InputParser(LH(i,:),Inputfile,MM,Minx,MSi,MFi);
        cmnd1=['"ANSYS\runwb2.exe" -F "RAF-TDMDO\EX1\MS\' Inputfile '.wbpj"'...
        ' -R "RAF-TDMDO\EX1\MS\' OutputFile '"'];
        tic
        system(cmnd1);
        A=xlsread('RAF-TDMDO\EX1\MS\Output.csv','A8:A100');
        f(i,j)=A(4);
        g(i,j)=norm([A(5) A(6) A(7)]);
        AA=textread(['RAF-TDMDO\EX1\MS\Test4_19_0\MD4_19_files\dp0\SYS-16\MECH\solve' '.out'],'%s');
        copyfile('RAF-TDMDO\EX1\MS\Test4_19_0\MD4_19_files\dp0\SYS-16\MECH\solve.out',...
            ['RAF-TDMDO\EX1\MS\out' num2str(h) '.out'])
        [r]=strcmp(AA,'LISTING');
        [ir]=find(r==1);
        [D]=find(r(ir+1:end,1)==1);
        K=[];
        disp(['Writing outputs for Run no. ' num2str(h) '...'])
        K = Write_Out(ir,AA);
        ys(i,j,:)=K;
        Es(i,j,:)=norm(ys(i,j,:)-ys0);
        inp(h,:)=[i LH(j,:)];
        h=h+1;
        
    end
end

for i=1:1:4
    for j=1:1:length(LH)
        Minx = [1 i]
        OutputFile=InputParser(LH(i,:),Inputfile,MM,Minx,MSi,MFi);
        cmnd1=['"ANSYS\runwb2.exe" -F "RAF-TDMDO\EX1\MF\' Inputfile '.wbpj"'...
        ' -R "RAF-TDMDO\EX1\MF\' OutputFile '"'];
        tic
        system(cmnd1);
        A=xlsread('RAF-TDMDO\EX1\MF\Output.csv','A8:G15');
        f(i,j)=A(4);
        g(i,j)=norm([A(5) A(6) A(7)]);
        if h>1
            yf(i,j,:)= A(8:15);
            Ef(i,j,:)=norm(yf(i,j,:)-yf0);
        end
        h=h+1;
    end
end



% close all
clear all
clc

%% Error space


ESS=fit(inp,Es);

ESF=fit(inp,Ef);

RAM_HAT0=[feval(ESS,[1,x])+feval(ESF,[1,x]) feval(ESS,[1,x])+feval(ESF,[2,x]) feval(ESS,[1,x])+feval(ESF,[3,x]) feval(ESS,[1,x])+feval(ESF,[4,x]);
    feval(ESS,[2,x])+feval(ESF,[1,x]) feval(ESS,[2,x])+feval(ESF,[2,x]) feval(ESS,[2,x])+feval(ESF,[3,x]) feval(ESS,[2,x])+feval(ESF,[4,x]);
    feval(ESS,[3,x])+feval(ESF,[1,x]) feval(ESS,[3,x])+feval(ESF,[2,x]) feval(ESS,[3,x])+feval(ESF,[3,x]) feval(ESS,[3,x])+feval(ESF,[4,x]);
    feval(ESS,[4,x])+feval(ESF,[1,x]) feval(ESS,[4,x])+feval(ESF,[2,x]) feval(ESS,[4,x])+feval(ESF,[3,x]) feval(ESS,[4,x])+feval(ESF,[4,x]);];

TR0=0.5;

lb=[0 0];
ub=[1 1];
bb_handle = @(x,MSi,MFi,Minx) ANSYSBB(x,MM,Minx);
[Xmin,fmin,output] = simple_mads2(X0,bb_handle,ESS,ESF,RAM_HAT0,MM,MSi,MFi,TR0,lb,ub)

end

function K = Write_Out(ir,AA)
if isempty(ir)
        disp('Check the Analysis for errors!')
else
    disp(['Running Run no. ' num2str(i) ' Finished Successfully!'])
    I=ir(1)+6;
    hh=I(1,1);
    Term=0;
%         hhh=1;
    for hhh=1:1:4      % index_n*D+I how to get n? : how many pages to read
%             for j=I(jj):6:48*2+I(jj)+2   % 48*2+jj+2 last line? : how many lines to read? 
            if strcmp((cell2mat(AA(hh,1))),'Set')==1
                break
                Term=1;
            end
            if strcmp((cell2mat(AA(hh,1))),'***')==1
                hh=hh+36;
            else
                K(hhh,1:2)=[str2num(cell2mat(AA(hh,1))) str2num(cell2mat(AA(hh+1,1)))];
%                     if hhh>1 && sum(k(hhh,:))==0
%                         K(hhh,:)=[];
%                     end
                hh=hh+2;
%                     hhh=hhh+1;
            end

%             end
    end
end
nz=find(K(:,1)==0);
K(nz,:)=[];
save(['KI_' num2str(i) '.mat'],'K')
end

function OutputFile=InputParser(x,InputFile,MM,Minx,MS,MF)

fid = fopen(InputFile);
OutputFile=['MI' num2str(index) '.wbpjn'];
data1 = textscan(fid, '%s', 'Delimiter', '\n', 'CollectOutput', true);
fclose(fid);
SearchString1=Var1;
SearchString2=Var2;
SearchString3=MSApp;
SearchString4=MSm;
SearchString5=MFm;
SearchString6=MFt;
SearchString7=MTC;

% SearchString8=Var8;
% modify the cell array
% find the position where changes need to be applied and insert new data
for I = 1:length(data1{1})
    tf1 = strcmp(data1{1}{I}, SearchString1); % search for this string in the array
    tf2 = strcmp(data1{1}{I}, SearchString2); % search for this string in the array
    tf3 = strcmp(data1{1}{I}, SearchString3); % search for this string in the array
    tf4 = strcmp(data1{1}{I}, SearchString4); % search for this string in the array
    tf5 = strcmp(data1{1}{I}, SearchString5); % search for this string in the array
    tf6 = strcmp(data1{1}{I}, SearchString6); % search for this string in the array
    tf7 = strcmp(data1{1}{I}, SearchString7); % search for this string in the array
    if tf1 == 1
        ReplaceString=num2str(x(1));
        data1{1}{I} = ReplaceString; % replace with this string
    end
    if tf2 == 1
        ReplaceString=num2str(x(2));
        data1{1}{I} = ReplaceString; % replace with this string
    end
    if tf3 == 1
        ReplaceString=MM.MS.APP{MS(Minx(1),2)};
        data1{1}{I} = ReplaceString; % replace with this string
    end
    if tf4 == 1
        ReplaceString=num2str(MM.MS.size{MS(Minx(1),3)});
        data1{1}{I} = ReplaceString; % replace with this string
    end
    if tf5 == 1
        ReplaceString=num2str(MM.MF.size{MF(Minx(2),2)});
        data1{1}{I} = ReplaceString; % replace with this string
    end
    if tf6 == 1
        ReplaceString=num2str(MM.MF.DT{MF(Minx(2),3)});
        data1{1}{I} = ReplaceString; % replace with this string
    end
    if tf7 == 1
        ReplaceString=num2str(MM.MF.DT{MF(Minx(2),3)});
        data1{1}{I} = ReplaceString; % replace with this string
    end
end

fid = fopen(OutputFile, 'w');
for I = 1:length(data1{1})
    fprintf(fid, '%s\n', char(data1{1}{I}));
end
fclose(fid);
end
