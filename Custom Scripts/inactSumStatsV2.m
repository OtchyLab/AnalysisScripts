function inactSumStatsV2(dataName, dataFolder,D)
%3/8/15
%Updated from the original version (inactSumStats.m) to remove interpolation of heat maps, reference the original experiment
%type, 
%Takes the output of the inactAnalysis scripts (the complete set of working data) and runs summary stats on the different
%conditions. All output data is saved to an on-disk mat-file that include bird names, dates, conditions, and results.

% clear all; close all;

%Set the output file location
outputLocation = 'C:\Users\Tim\Desktop\Nif Project Figures\Inactivations\';
outName = 'Inactivation Summary Stats EMD 05142015.mat';

%Select the dataset to use
%[dataName, dataFolder] = uigetfile('C:\Users\Tim\Desktop\Nif Project Figures\Inactivations\Grn114\*.mat');
% dataFolder = 'C:\Users\Tim\Desktop\Nif Project Figures\Inactivations\Grn108';
% dataName = 'Grn108_140828_inactDataset.mat';
load([dataFolder, filesep, dataName]);
close all; %close the windows

%Sets the manipulation type based on the original dataset values
injType = expType;

%Parse name to grab bird name, date and type of experiment
sp = regexp(dataName,'_', 'split');
output.name = char(sp{1});
output.date = char(sp{2});
output.filename = [dataFolder, dataName];
output.type = injType;

%Comparison labels
output.labels = [{'pre-inj'}, {'inj-post'}, {'pre-post'}];

%%%%%%%%%%%%%%%%%%%%%
% Capture distributions for later plotting
%%%%%%%%%%%%%%%%%%%%%
%1-D
output.preY = preY;
output.injY = injY;
output.postY = postY;

%2-D
%Capture originals
% output.preNcnts = preNcnts;
% output.injNcnts = injNcnts;
% output.postNcnts = postNcnts;

output.preN = preN;
output.injN = injN;
output.postN = postN;

%Filter the original PDFs with a gaussian and store in original variable for ease
output.filter.size = 7;
output.filter.sigma = 3;
output.filter.h = fspecial('gaussian', output.filter.size, output.filter.sigma);

output.preNsm = imfilter(preN,output.filter.h,'replicate');
output.injNsm = imfilter(injN,output.filter.h,'replicate');
output.postNsm = imfilter(postN,output.filter.h,'replicate');

%Copy out for ease of coding
preN = output.preNsm;
injN = output.injNsm;
postN = output.postNsm;

thetaPI = calcEMD(preN,injN,D); %Swapped for EMD
thetaIP = calcEMD(injN,postN,D);
thetaPP = calcEMD(preN,postN,D);
output.D2emd = [thetaPI, thetaIP, thetaPP];


%%%%%%%%%%%%%%%%%%%%%
% Save to output
%%%%%%%%%%%%%%%%%%%%%
m = exist([outputLocation, outName]);
if m ==2 
    %File already exists
    clear('inactStats');
    load([outputLocation, outName], 'inactStats')
    inactStats(end+1) = output;
else
    %No file yet created
    inactStats = output;
end

% Save the updated data to file
save([outputLocation, outName], 'inactStats')

clear all

display('done')













