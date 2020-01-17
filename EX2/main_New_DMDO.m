%% Clean directory 
% clear all
clc
PATH.temp='Temp';
dinfo = dir(PATH.temp);
dinfo([dinfo.isdir]) = [];   %skip directories
filenames = fullfile(PATH.temp, {dinfo.name});
% if ~cellfun(@isempty, filenames)
%     delete(filenames{:})
% end
PATH.Working=pwd;

savecache=1;
%% Define the DFD for initial analysis
i=0;
PATH.FluidRes='D:\Data\Research\Publications\SMO\ANSYS files\MMM2\MMM01_files\dp0\CFX\CFX';
PATH.Cache='D:\Data\Research\Publications\SMO\ANSYS files\MMM2\cache\dp';
PATH.Proj='D:\Data\Research\Publications\SMO\ANSYS files\MMM2';
PATH.MECHRes='D:\Data\Research\Publications\SMO\ANSYS files\MMM2\MMM01_files\dp0\SYS\MECH';
if savecache==1
    mkdir([PATH.Cache num2str(i)]);
    copyfile([PATH.FluidRes],[PATH.Cache num2str(i)]);
end
SearchString=[PATH.FluidRes 'Fluid Flow CFX_001.res'];

for i=10:1:10
    copyfile([PATH.Working '\Trial_' num2str(0) '.wbjn'],...
        [PATH.temp '\Trial_' num2str(i-1) '.wbjn'])
    InputFile=[PATH.temp '\Trial_' num2str(i-1) '.wbjn'];
    OutputFile=['Trial_' num2str(i-1) '.wbjn'];
    fid = fopen(InputFile);
    data1 = textscan(fid, '%s', 'Delimiter', '\n', 'CollectOutput', true);
    fclose(fid);
    % modify the cell array
    % find the position where changes need to be applied and insert new data
    for I = 1:length(data1{1})
        tf = strcmp(data1{1}{I}, SearchString); % search for this string in the array
        if tf == 1
            ReplaceString=[PATH.FluidRes '\Fluid Flow CFX_00' num2str(1) '.res'];
            data1{1}{I} = ReplaceString; % replace with this string
        end
    end
    % write the modified cell array into the text file
    fid = fopen(OutputFile, 'w');
    for I = 1:length(data1{1})
        fprintf(fid, '%s\n', char(data1{1}{I}));
    end
    fclose(fid);
    disp(['Running Run no. ' num2str(i) '...'])
    cmnd1=['"runwb2.exe" -F "' PATH.Proj '\MMM01.wbpj"'...
        ' -R "' PATH.temp '\' OutputFile '"'];
    tic
    [status1,cmdout1] = system(cmnd1,'-echo')
    disp(['Postprocess Run no. ' num2str(i) '...'])
    InputFile2=[PATH.Proj '\MMM01_files\dp0\SYS\MECH\ds' '.dat'];
    copyfile([PATH.Proj '\MMM01_files\dp0\SYS\MECH\ds' '.dat'],...
        [PATH.temp 'ds' num2str(i) '.dat'])
    fid22 = fopen(InputFile2);
    data2 = textscan(fid22, '%s', 'Delimiter', '\n', 'CollectOutput', true);
    
    fclose(fid22);
    Nodes=zeros(length(data2{1}),4);
    SearchString=['NSOL,2,']; 
    for I = 1:length(data2{1})
        if ~isempty(str2num(data2{1}{I})) && size(str2num(data2{1}{I}),2)==4
            Nodes(I,:) = (str2num(data2{1}{I})); % search for this string in the array
        end
    end
    
    nr0=find(Nodes(:,1)==1);
    Nodes(nr0,:)=[];
    [nr1,nc1]=find(Nodes(:,4)<-1);
    [nr2,nc2]=find(Nodes(nr1,2)==max(Nodes(nr1,2)));
    [nr3,nc3]=find(Nodes(nr1(nr2),3)==max(abs(Nodes(nr1(nr2),3))));
%     [nr,nc]=find(Nodes(nr1,3)==max(Nodes(nr1,3)));
    ND(i,1)=Nodes(nr1(nr2(nr3(1,1))),1);
%     tff=0;
%     for I = 1:length(data2{1})
%         if size(data2{1}{I},2)>=7
%             tff = strcmp(data2{1}{I}(1:7), SearchString); % search for this string in the array
%         end
%         if tff == 1
%             ReplaceString=['NSOL,2,' num2str(ND) ',U,Y, UY_2, '];
%             data2{1}{I} = ReplaceString; % replace with this string
%         end
%     end
%     fid22 = fopen(InputFile2, 'w');
%     for I = 1:length(data2{1})
%         fprintf(fid22, '%s\n', char(data2{1}{I}));
%     end
%     fclose(fid);
    % write the modified cell array into the text file
    InputFile3=['D:\Data\post' '.txt'];
%     OutputFile=['Trial_' num2str(i) '.wbjn'];
    fid = fopen(InputFile3);
    data3 = textscan(fid, '%s', 'Delimiter', '\n', 'CollectOutput', true);
    fclose(fid);
    % modify the cell array
    % find the position where changes need to be applied and insert new data
    ReplaceString=['NSOL,2,' num2str(ND(end,1)) ',U,Y, UY_2, '];
    data3{1}{30} = ReplaceString; % replace with this string
            
    % write the modified cell array into the text file
    fid = fopen(InputFile3, 'w');
    for I = 1:length(data3{1})
        fprintf(fid, '%s\n', char(data3{1}{I}));
    end
    fclose(fid);
    
    
    
% % % % % %     cmnd2=['"runwb2.exe" -F "' PATH.Proj '\MMM01.wbpj"'...
% % % % % %         ' -R "' PATH.Working '\TC2.wbjn"']
% % % % % %     
% % % % % %     
% % % % % %     [status2,cmdout2] = system(cmnd2,'-echo');
% % % % % %     ICC(i,1)=toc;
% % % % % %     disp(['Time Elapsed for Run no. ' num2str(i) '...' num2str(ICC(i,1))])
    
    
    AA=textread([PATH.Proj '\MMM01_files\dp0\SYS\MECH\solve' '.out'],'%s');
    copyfile([PATH.MECHRes '\solve.out'],...
        [PATH.temp '\' num2str(i) '.out'])
    [r]=strcmp(AA,'LISTING');
    [ir]=find(r==1);
    [D]=find(r(ir+1:end,1)==1);
    K=[];
    disp(['Writing outputs for Run no. ' num2str(i) '...'])

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
    if sum(K(:,2))<=1E-9
        disp('Check the Analysis for errors!')
%         break
    end
    save([PATH.temp '\KI_' num2str(i) '.mat'],'K')
end


