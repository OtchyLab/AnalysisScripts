% function electroSumStatsV2
%4/17/15
%Updated from the original version (electroSumStats.m) to remove interpolation of heat maps, reference
%the original experiment type. Principally for the ELECTROLYTIC Nif lesions.
%Takes the output of the electroAnalysis.m scripts (across several days, and all in the same folder) 
%and runs summary stats on the different conditions. All output data is saved to an on-disk mat-file that include
%bird names, dates, conditions, and results.

%Set the output file location
outputLocation = 'C:\Users\Tim\Desktop\Nif Project Figures\ElectroLesions\060315\';
outName = 'Electro Summary Stats days3-4.mat';

%Select the dataset to use
%[dataName, dataFolder] = uigetfile('C:\Users\Tim\Desktop\Nif Project Figures\Inactivations\Grn114\*.mat');
% dataFolder = 'C:\Users\Tim\Desktop\Nif Project Figures\Pharma\Grn089\';
% dataName = 'Grn089_140625_pharmaDataset.mat';
% load([dataFolder, filesep, dataName]);
% close all; %close the windows

%Ask user for data directory
dataFolder = uigetdir('C:\Users\Tim\Desktop\Nif Project Figures\ElectroLesions\060315\');

%Get all data files that match the criterion
crit = '*electroDataset.mat';
filelist = dir([dataFolder, filesep, crit]);

%Cycle through the filelist and sequentially process the datasets from the whole folder
data = [];
output = [];
for i = 1:length(filelist)
    %Select the file name
    dataName = filelist(i).name;

    %Load the data file
    t = load([dataFolder, filesep, dataName]);
    close all; %close the windows
    
    %Parse name to grab bird name, date and type of experiment
    sp = regexp(dataName,'_', 'split');
    output.name = char(sp{1});
    idx = t.offset+5; %add 5 to avoid the zero/negative offsets; 5 = lesion day.
    
    %Correct for L-5d to L-4d switch so rest of code works
%     if strcmp(char(sp{1}), 'Grn091') && t.conIdx==2
%         t.conIdx = 1;
%         t.offset = -5;
%         display('Made the Grn091 switch!')
%     end
    
    output.date{idx} = char(sp{2});
    output.filename{idx} = [dataFolder, filesep, dataName];
%     output.type{idx} = t.conditions{t.conIdx};
    output.offset(idx) = t.offset;
%     output.conditions = t.conditions;
    
    % Capture distributions for later stats
    output.mornY(idx,:) = t.mornY;
    output.nightY(idx,:) = t.nightY;
    output.allY(idx,:) = t.allY;
    
    %Filter the original PDFs with a gaussian and store in original variable for ease
    output.filter.size = 7;
    output.filter.sigma = 3;
    output.filter.h = fspecial('gaussian', output.filter.size, output.filter.sigma);
    
    output.mornN(idx,:,:) = imfilter(t.mornN,output.filter.h,'replicate');
    output.nightN(idx,:,:) = imfilter(t.nightN,output.filter.h,'replicate');
    output.allN(idx,:,:) =imfilter(t.allN,output.filter.h,'replicate');
    
    if isfield(t, 'postY') && isfield(t, 'postN')
        output.postY(idx,:) = imfilter(t.postY,output.filter.h,'replicate');
        output.postN(idx,:,:) = imfilter(t.postN,output.filter.h,'replicate');
    end
    
%     output.mornN(idx,:,:) = t.mornN;
%     output.nightN(idx,:,:) = t.nightN;
%     output.allN(idx,:,:) = t.allN;
%     
%     if isfield(t, 'postY') && isfield(t, 'postN')
%         output.postY(idx,:) = t.postY;
%         output.postN(idx,:,:) = t.postN;
%     end
    

    %Clear the dataset
    clear('t')
    
end

%%%%%%%%%%%%%%%
% Averaged Pre
%%%%%%%%%%%%%%%
%Hard coded matrix of all the comarisons to be made
cntPntr = find(output.offset~=0,1,'first');
lastPntr = find(output.offset~=0,1,'last');

% cmpMat = [(output.mornN(cntPntr,:,:));...         %control all day
%                  (output.postN(5,:,:));...               %post lesion
%                  (output.nightN(5,:,:));...              %night of lesion
%                  (output.mornN(6,:,:));...                  %1 day after lesion
%                  (output.mornN(7,:,:));...                  %2 day after lesion
%                  (output.mornN(lastPntr,:,:))];           %last day lesion

cmpMat = [(output.mornN(8,:,:));...                  %3 day after lesion
                 (output.mornN(9,:,:))];           %4 day after lesion

%             check this matrix before running
             
%Copy to output structure
output.cmpMat = cmpMat;

%Load the cost matrix
load('CostMatrix 281x161.mat');

%Loop through the comparisons and calculate stats
for i = 1:size(cmpMat,1)
    
    if isempty(squeeze(cmpMat(i, :, :)))
        %If this subset is empty, do nothing???
        
    else
        %Otherwise, do the processing on thie dataset.

        %%%%%%%%%%%%%%%%%%%%
        % 2-D measures
        %%%%%%%%%%%%%%%%%%%%
        %Simplify variables for reading
        N1 = squeeze(output.mornN(5,:,:)); % The baseline distribution
        output.baseN = N1;

           %Simplify variables for reading
        N2 = squeeze(cmpMat(i, :, :));
        i
        output.D2emd(i) = calcEMD(N1, N2, D);
    end
end

%%%%%%%%%%%%%%%%%%%%%
% Save to output
%%%%%%%%%%%%%%%%%%%%%
 m = exist([outputLocation, outName]);
if m ==2 
    %File already exists
    clear('electroStats');
    load([outputLocation, outName], 'electroStats')
    electroStats(end+1) = output;
else
    %No file yet created
    electroStats = output;
end

%Save the updated data to file
save([outputLocation, outName], 'electroStats')

display(['done with ' output.name])

clear all



