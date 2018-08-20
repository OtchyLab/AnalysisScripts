function processPartialFolder

%Ask user for data directory
dataFolder = uigetdir('C:\Users\Tim\Desktop\Nif Project Figures\ElectroLesions\Lesion Size-Effect\');

%Get all data files that match the criterion
crit = '*allData.mat';
filelist = dir([dataFolder, filesep, crit]);

%Cycle through the filelist and sequentially process the datasets
for i = 1:length(filelist)
    dataName = filelist(i).name;
    
    %Generate summary stats
    getoutput(dataName, dataFolder)
end

display(['Finished with folder: ' dataFolder])



function getoutput(dataName, dataFolder)
%Load data from file
load([dataFolder filesep dataName])
close all

output = [];
output.name = sp{1};
output.date = sp{2};
output.filename = annotName;
output.chan = str2num(sp{3}(end-4));
output.lesionFrac = lesionFrac;

    %HVC mean activity
    m = []; s = [];
m(1) = mean(meanPow(preInx));
m(2) = mean(meanPow(postInx));
m(3) = mean(meanPow(nightInx));
s(1) = std(meanPow(preInx));
s(2) = std(meanPow(postInx));
s(3) = std(meanPow(nightInx));
output.powMeans = m; output.powStd = s;

    %HVC recovery
    m = []; s = [];
m(1) = mean(recovCorr(preInx));
m(2) = mean(recovCorr(postInx));
m(3) = mean(recovCorr(nightInx));
s(1) = std(recovCorr(preInx));
s(2) = std(recovCorr(postInx));
s(3) = std(recovCorr(nightInx));
output.recoveryMeans = m; output.recoveryStd = s;

%%%%%%%%%%%%%%%%%%%%%
% Save to output location
%%%%%%%%%%%%%%%%%%%%%

outputLocation = 'C:\Users\Tim\Desktop\Nif Project Figures\ElectroLesions\';
outPerm = 'Partial Lesion Summary Data XXX.mat';
m = exist([outputLocation, outPerm]);
if m ==2 
    %File already exists
    clear('partialSumStats');
    load([outputLocation, outPerm], 'partialSumStats')
    partialSumStats(end+1) = output;
else
    %No file yet created
    partialSumStats = output;
end

%Save the updated data to file
save([outputLocation, outPerm], 'partialSumStats')
