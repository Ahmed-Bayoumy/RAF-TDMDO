function [bbo] = ANSYSBB(Xtry,MM,Minx,MS,MF)
%ANSYSBB Summary of this function goes here
%   Detailed explanation goes here
OutputFile=InputParser(Xtry,'MI.wbjn',MM,Minx,MS,MF);
cmnd1=['"ANSYS\runwb2.exe" -F "RAF-TDMDO\EX1\MS\' Inputfile '.wbpj"'...
' -R "RAF-TDMDO\EX1\MS\' OutputFile '"'];
tic
[status1,cmdout1] = system(cmnd1);
A=xlsread('RAF-TDMDO\EX1\MS\Output.csv','A8:A100');
AA = xlsread('RAF-TDMDO\EX1\MS\Output.csv','A8:G8');
f=AA(5);
g=norm([A(4)-3E6 A(6)-0.3 A(7)]);

bbo=[f g];



end



function OutputFile=InputParser(x,InputFile,MM,Minx,MS,MF)

fid = fopen(InputFile);
OutputFile=[InputFile num2str(index)];
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

% function K = Write_Out(ir,AA)
% if isempty(ir)
%         disp('Check the Analysis for errors!')
% else
%     disp(['Running Run no. ' num2str(i) ' Finished Successfully!'])
%     I=ir(1)+6;
%     hh=I(1,1);
%     Term=0;
% %         hhh=1;
%     for hhh=1:1:4      % index_n*D+I how to get n? : how many pages to read
% %             for j=I(jj):6:48*2+I(jj)+2   % 48*2+jj+2 last line? : how many lines to read? 
%             if strcmp((cell2mat(AA(hh,1))),'Set')==1
%                 break
%                 Term=1;
%             end
%             if strcmp((cell2mat(AA(hh,1))),'***')==1
%                 hh=hh+36;
%             else
%                 K(hhh,1:2)=[str2num(cell2mat(AA(hh,1))) str2num(cell2mat(AA(hh+1,1)))];
% %                     if hhh>1 && sum(k(hhh,:))==0
% %                         K(hhh,:)=[];
% %                     end
%                 hh=hh+2;
% %                     hhh=hhh+1;
%             end
% 
% %             end
%     end
% end
% nz=find(K(:,1)==0);
% K(nz,:)=[];
% save(['KI_' num2str(i) '.mat'],'K')
% end
