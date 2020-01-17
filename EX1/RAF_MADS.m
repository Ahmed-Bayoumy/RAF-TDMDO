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
%  https://github.com/ahmedbayoumy/RAF_MADS                                           %
%-------------------------------------------------------------------------------------%


function [Xmin,fmin,output] = RAF_MADS(X0,bb_handle,ESS,ESF,RAM_HAT0,MM,MSi,MFi,TR0,lb,ub)

% Default options
default_options.budget = +inf;
default_options.tol    = 1e-9;
default_options.psize_init = 1.0;
default_options.display         = false;
default_options.opportunistic   = true;
default_options.check_cache     = true;
default_options.Inputfile     = 'MI.wbpjn';
Inputfile = default_options.Inputfile ; 
% Options
if ~exist('options','var')
    options = struct;
end
names = fieldnames(default_options);
for i = 1:length(names)
    fieldstr = names{i};
    if ~isfield(options,fieldstr)
        options.(fieldstr) = default_options.(fieldstr);
    end
end

DIM = size(X0,2);
if length(lb)~=DIM
    error('Wrong lb dimension');
end
if length(ub)~=DIM
    error('Wrong ub dimension');
end
bbe = 0;
CACHE = [];
fCACHE=[];
fmin = +inf;
hmin = +inf;

Rep=[];
% Eval=0;
% Evaluation of starting point(s)
    




fail=0;
% Variable scaling
scaling = (ub-lb)/10.0;
scaling(isinf(scaling)) = 1.0;
scaling = diag(scaling);

% Init
dir_success = [];
psize = options.psize_init;
psize_success = psize;
psize_max = 0;
iter = 1;
success=false;
TRF = 0;
TRk=TR0;
while true
    %% Search step (optional)
    if TRF==1
        N=2;
        px0 = 6;
        X0 = lhsdesign(px0,N);
        for i=1:size(X0,1)
                Xtry = X0(i,:);
                
                RAM_HATk=[feval(ESS,[1,Xtry])+feval(ESF,[1,Xtry]) feval(ESS,[1,Xtry])+feval(ESF,[2,Xtry]) feval(ESS,[1,Xtry])+feval(ESF,[3,Xtry]) feval(ESS,[1,Xtry])+feval(ESF,[4,Xtry]);
            feval(ESS,[2,Xtry])+feval(ESF,[1,Xtry]) feval(ESS,[2,Xtry])+feval(ESF,[2,Xtry]) feval(ESS,[2,Xtry])+feval(ESF,[3,Xtry]) feval(ESS,[2,Xtry])+feval(ESF,[4,Xtry]);
            feval(ESS,[3,Xtry])+feval(ESF,[1,Xtry]) feval(ESS,[3,Xtry])+feval(ESF,[2,Xtry]) feval(ESS,[3,Xtry])+feval(ESF,[3,Xtry]) feval(ESS,[3,Xtry])+feval(ESF,[4,Xtry]);
            feval(ESS,[4,Xtry])+feval(ESF,[1,Xtry]) feval(ESS,[4,Xtry])+feval(ESF,[2,Xtry]) feval(ESS,[4,Xtry])+feval(ESF,[3,Xtry]) feval(ESS,[4,Xtry])+feval(ESF,[4,Xtry]);];
                LAMBDA = (RAM_HATk-RAM_HAT0);
                LAMBDAn=[norm(LAMBDA(1,1)) norm(LAMBDA(1,2)) norm(LAMBDA(1,3)) norm(LAMBDA(1,4));
                    norm(LAMBDA(2,1)) norm(LAMBDA(2,2)) norm(LAMBDA(2,3)) norm(LAMBDA(2,4));
                    norm(LAMBDA(3,1)) norm(LAMBDA(3,2)) norm(LAMBDA(3,3)) norm(LAMBDA(3,4));
                    norm(LAMBDA(4,1)) norm(LAMBDA(4,2)) norm(LAMBDA(4,3)) norm(LAMBDA(4,4))];

                [r,c]=find(LAMBDAn<=TRk)
                
                r=((ffk+hhk)-(ff+hh))/((fk+hk)-(f+h));

                if sum(any(LAMBDAn<=TRk))
                    ES=feval(ESS,[r,Xtry])
                    EF=feval(ESF,[c,Xtry])
                    Minx = [r c];
                    
                else
                    for i=1:1:size(X0,1)
                        XXtry = X0(i,:);
                        [ES,EF,RAM_HATk] = Re_evaluate(Inputfile,MM,MSi,MFi);
                    end
                end
                
                

                %disp(num2str(Xtry))

                if any(Xtry>ub) || any(Xtry<lb)
                    Xtry
                    ub
                    lb
                    error('Starting point is not within the box constraints');
                end

                % Evaluation
                bbe = bbe+1;
                bbo = bb_handle(Xtry,MM,Minx,MSi,MFi);

                % eval_fh
                while isempty(bbo)
                    Xtry=Xtry+((Xtry-lb)/10);
                    X0(i,:)=Xtry;
                    bbo = bb_handle(Xtry,MM,Minx,MSi,MFi);
        %             initial=initial-1;
                end
                [ftry,htry] = eval_fh(bbo);
                
                

                % Add to the CACHE
                if options.check_cache
        %             fCACHE
                    CACHE(end+1,:) = Xtry;
                    fCACHE(end+1,:)=[ftry htry];
        %             fCACHE
                end

                % Improvement of the feasibility or of the objective ?
                success = ( (hmin>0) && (htry<hmin) ) || ( (htry==0) && (ftry<fmin) );
                r=((ffmin+hhmin)-(fftry+hhtry))/((fmin+hmin)-(ftry+htry));
                if r<R1
                    TRk=c1*TRk;
                elseif r>R2
                    TRk=c2*TRk;
                else

                end
                % If success, save values
                if success
                    fmin = ftry;
                    hmin = htry;
                    Xmin = Xtry;
                end
                if options.display
                    disp(['x0     : ' num2str(ftry,4) ' (hmin = ' num2str(htry,4) ')']);
                end

                RAM_HAT0=RAM_HATk;

        end
    else
    
    %% Poll step
%     Eval=0;
    % Build polling directions
        msize = min(psize^2,psize);
        rho = psize/msize;
        % Generate direction
        if success
            v=dir_success;
        else
            v = randn(DIM,1);
        end
        % Normalize
        v = v/norm(v);
        % Build Householder matrix
        H = eye(DIM)-2*v*v';
        % Normalization of each column
        H = H*diag(max(abs(H)).^-1);
        % Rounding (and transpose)
        H = msize*ceil(rho*H)';
        % Take the opposite directions
        H = [H;-H];

        % scaling of the directions
        H = H*scaling;

        % Build POLL / central point

        POLL = ones(2*DIM,1)*Xmin + H;

        % Shuffle poll
        POLL = POLL(randperm(2*DIM),:);

        % Evaluate points of the poll
        success = false;
        Xold = Xmin;
        for i=1:size(POLL,1)

            % Get point
            Xtry = POLL(i,:);

            % Snap to bounds
            Xtry = max(Xtry,lb);
            Xtry = min(Xtry,ub);

            % Search in CACHE
            if options.check_cache && ~isempty(CACHE)
                DC = CACHE-ones(size(CACHE,1),1)*Xtry;
                DC = abs(DC)<options.tol;
                if any(all(DC,2))
                    continue;
                end
            end

            % Evaluation
            bbe = bbe+1;
            bbo = bb_handle(Xtry,MM,Minx,MSi,MFi);



            % compute h and f
            [ftry,htry] = eval_fh(bbo);

            if options.check_cache && ~isempty(fCACHE)
                DC = CACHE-ones(size(CACHE,1),1)*Xtry;
    %             isinf
                [Val,Indx]=min((abs(DC(:,1))+abs(DC(:,2))));
                FDC1 = fCACHE(Indx,:)-[ftry htry];
                FDC = ((abs(FDC1)./abs(fCACHE(Indx,:)))>0.3);

    %             Val
                if ~isinf(FDC1(1,1)) && any((FDC)) && Val(1,1)<=0.05
                    Rep(end+1,:)=[iter Indx Xtry fCACHE(Indx,:) ftry htry fCACHE(Indx,:)-[ftry htry] (abs(FDC)./abs(fCACHE(Indx,:))) FDC];
                    save('Rep.mat','Rep')
    %                 Eval=1;
                    psize=psize*1.1;
                    continue;
                end
            end

            % Add to the CACHE
            if options.check_cache
                fCACHE
                CACHE(end+1,:) = Xtry
                fCACHE(end+1,:)=[ftry htry];
            end

            % Save values
            if ( (htry==0) && (ftry<fmin) )
                success = true;
                r=((fftry+hhtry)-(ffmin+hhmin))/((fmin+fmin)-(ftry+htry));
                if r<R1
                    TRk=c1*TRk;
                elseif r>R2
                    TRk=c2*TRk;
                end
                
                fmin = ftry;
                hmin = htry;
                
                if options.display
                    disp(['Succes : ' num2str(fmin) ' (hmin = ' num2str(hmin) ')']);
                end
                Xmin = Xtry;
    %             dist=Xmin-Xreg
                dir_success = Xtry-Xold;
                psize_success = psize;
                psize_max = max(psize,psize_max);
            end

            if bbe>=options.budget
                break; % Reached the total budget
            end
            if success && options.opportunistic
                break; % Quit the evaluation of the POLL
            end

        end

        if success
            psize = psize*2;
    %         x_last_success=Xmin;
            fail=0;
        else
            psize = psize/2;
            fail=fail+1;
        end

    %     if fail>4
    %         Xtry=Xmin;
    %     end

        if options.display
            disp(['iter=' num2str(iter) ' bbe=' num2str(bbe) ' psize=' num2str(psize,3) '  hmin=' num2str(hmin,3) '  fmin=' num2str(fmin,3)  '  htry=' num2str(htry,3)]);
            Xtry
        end

        if (abs(psize)<options.tol) || (bbe>=options.budget)
            break;
        end
        
        
        iter = iter+1;
        RAM_HAT0=RAM_HATk;
    end
end


if options.display
    disp(['mads break - iter ' num2str(iter,2) ', psize ' num2str(psize,2) ', bbe ' num2str(bbe)]);
    disp(['Final  : ' num2str(fmin) ' (hmin = ' num2str(hmin) ')']);
end

% Build the output
output.fmin = fmin;
output.hmin = hmin;
output.psize = psize;
output.psize_success = psize_success;
output.psize_max = psize_max;

end % end function simple_mads



function [f,h] = eval_fh(bbo)
  % Concatenation of bbo from f & c functions
  f = bbo(1);
  if length(bbo)==1
      h = 0;
  else
      c = bbo(2:end);
      % Check non valid values
      if any(isnan(c))
          h = +inf;
      else
          h = sum(max(c,0).^2);
      end
  end
  % Penalize objective
  if isnan(f) || (h>0)
    f = +inf;
  end
end % end function eval_fh

function [f,g,ESS,ESF,RAM_HAT] = Re_evaluate(Inputfile,MM,MSi,MFi)
for i=1:1:4
    for j=1:1:length(XXtry)
        Minx = [i 2]
        OutputFile=InputParser(XXtry(i,:),Inputfile,MM,Minx,MSi,MFi);
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
        inp(h,:)=[i XXtry(j,:)];
        h=h+1;
        
    end
end

for i=1:1:4
    for j=1:1:length(XXtry)
        OutputFile=InputParser(XXtry(i,:),Inputfile,index);
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

