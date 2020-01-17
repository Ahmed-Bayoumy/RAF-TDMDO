%% Clean directory 
PATH.temp='Temp';
dinfo = dir(PATH.temp);
dinfo([dinfo.isdir]) = [];   %skip directories
filenames = fullfile(PATH.temp, {dinfo.name});
if ~cellfun(@isempty, filenames)
    delete(filenames{:})
end
clear all
clc

%% Define the Analysis 

SearchString='File Name = C:\Data\Project files\Flutter\NX\test MD\Test4_19_0\MD4_19_files\dp0\CFX-1\CFX\Fluid Flow CFX_003.res';

for i=187:1:300
    copyfile('C:\software\sumo-toolbox\Prelimres\TooSimple\bspline\SQP\optimize\slp_sqp\MFM\slp_sqp\Ch07\code\Newfolder\MFM\Trial_C.wbjn',...
        ['C:\software\sumo-toolbox\Prelimres\TooSimple\bspline\SQP\optimize\slp_sqp\MFM\slp_sqp\Ch07\code\Newfolder\MFM\Trial_' num2str(i) '.wbjn'])
    InputFile=['C:\software\sumo-toolbox\Prelimres\TooSimple\bspline\SQP\optimize\slp_sqp\MFM\slp_sqp\Ch07\code\Newfolder\MFM\Trial_' num2str(i) '.wbjn'];
    OutputFile=['Trial_' num2str(i) '.wbjn'];
    fid = fopen(InputFile);
    data1 = textscan(fid, '%s', 'Delimiter', '\n', 'CollectOutput', true);
    fclose(fid);
    % modify the cell array
    % find the position where changes need to be applied and insert new data
    for I = 1:length(data1{1})
        tf = strcmp(data1{1}{I}, SearchString); % search for this string in the array
        if tf == 1
            ReplaceString=['File Name = C:\Data\Project files\Flutter\NX\test MD\Test4_19_0\MD4_19_files\dp0\CFX-1\CFX\Fluid Flow CFX_' num2str(i+68) '.res'];
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
    cmnd1=['"C:\Program Files\ANSYS Inc\v190\Framework\bin\Win64\runwb2.exe" -B -F "C:\Data\Project files\Flutter\NX\test MD\Test4_19_0\MD4_19.wbpj"'...
        ' -R "C:\software\sumo-toolbox\Prelimres\TooSimple\bspline\SQP\optimize\slp_sqp\MFM\slp_sqp\Ch07\code\Newfolder\MFM\' OutputFile '"'];
    tic
    [status1,cmdout1] = system(cmnd1);
    disp(['Postprocess Run no. ' num2str(i) '...'])
    InputFile2=['C:\Data\Project files\Flutter\NX\test MD\Test4_19_0\MD4_19_files\dp0\SYS-16\MECH\ds' '.dat'];
    copyfile('C:\Data\Project files\Flutter\NX\test MD\Test4_19_0\MD4_19_files\dp0\SYS-16\MECH\ds.dat',...
        ['C:\software\sumo-toolbox\Prelimres\TooSimple\bspline\SQP\optimize\slp_sqp\MFM\slp_sqp\Ch07\code\Newfolder\MFM\DS\ds' num2str(i) '.dat'])
    fid = fopen(InputFile2);
    data2 = textscan(fid, '%s', 'Delimiter', '\n', 'CollectOutput', true);
    fclose(fid);
    Nodes=zeros(length(data2{1}),4);
    for I = 1:length(data2{1})
        if ~isempty(str2num(data2{1}{I})) && size(str2num(data2{1}{I}),2)==4
            Nodes(I,:) = (str2num(data2{1}{I})); % search for this string in the array
        end
    end
    
    nr0=find(Nodes(:,2)>15);
    Nodes(nr0,:)=[];
    [nr1,nc1]=find(Nodes(:,2)==max(Nodes(:,2)));
%     [nr,nc]=find(Nodes(nr1,3)==max(Nodes(nr1,3)));
    ND=Nodes(nr1(1,1),1);
    % write the modified cell array into the text file
    InputFile3=['C:\Data\post' '.txt'];
%     OutputFile=['Trial_' num2str(i) '.wbjn'];
    fid = fopen(InputFile3);
    data3 = textscan(fid, '%s', 'Delimiter', '\n', 'CollectOutput', true);
    fclose(fid);
    % modify the cell array
    % find the position where changes need to be applied and insert new data
    ReplaceString=['NSOL,2,' num2str(ND) ',U,Y, UY_2, '];
    data3{1}{30} = ReplaceString; % replace with this string
            
    % write the modified cell array into the text file
    fid = fopen(InputFile3, 'w');
    for I = 1:length(data3{1})
        fprintf(fid, '%s\n', char(data3{1}{I}));
    end
    fclose(fid);
    
    
    
    cmnd2=['"C:\Program Files\ANSYS Inc\v190\Framework\bin\Win64\runwb2.exe" -B -F "C:\Data\Project files\Flutter\NX\test MD\Test4_19_0\MD4_19.wbpj"'...
        ' -R "C:\software\sumo-toolbox\Prelimres\TooSimple\bspline\SQP\optimize\slp_sqp\MFM\slp_sqp\Ch07\code\Newfolder\MFM\TC2.wbjn"'];
    
    
    [status2,cmdout2] = system(cmnd2);
    ICC(i,1)=toc;
    disp(['Time Elapsed for Run no. ' num2str(i) '...' num2str(ICC(i,1))])
    
    
    AA=textread(['C:\Data\Project files\Flutter\NX\test MD\Test4_19_0\MD4_19_files\dp0\SYS-16\MECH\solve' '.out'],'%s');
    copyfile('C:\Data\Project files\Flutter\NX\test MD\Test4_19_0\MD4_19_files\dp0\SYS-16\MECH\solve.out',...
        ['C:\software\sumo-toolbox\Prelimres\TooSimple\bspline\SQP\optimize\slp_sqp\MFM\slp_sqp\Ch07\code\Newfolder\MFM\Solve\solve' num2str(i) '.out'])
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
    save(['KI_' num2str(i) '.mat'],'K')
end


