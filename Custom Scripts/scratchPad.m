%Script to calculate the change in acoustic power over recovery
clear all

%Set the output file location
outLocation = 'C:\Users\Tim\Desktop\Nif Project Figures\Song Power Recovery\';
outName = 'Power Recovery XXX.mat';

%Ask user for data directory
dataFolder = uigetdir('C:\Users\Tim\Desktop\Electro Songs\');

%Get all data files that match the criterion
crit = '*Songs.mat';
filelist = dir([dataFolder, filesep, crit]);

%Cycle through the filelist and sequentially parse the datasets from the whole folder
runPower = [];
runTime = [];
for i = 1:length(filelist)
    %Select the file name
    dataName = filelist(i).name;
    
    %Load the data file
    motifSpecs = []; filenames = [];
    load([dataFolder, filesep, dataName], 'filt_audio', 'filenames');
    
    %Parse name to grab bird name, date
    sp = regexp(dataName,' ', 'split');
    birdname = char(sp{1});
    dates = char(sp{2});
    
    if strcmp(birdname, 'Grn046')
        critTime = datenum(2014, 5, 15, 11, 32, 00); %5/15/2014 @ 11:32
    elseif strcmp(birdname, 'Grn121')
        critTime = datenum(2014, 9, 6, 14, 5, 00); %9/6/2014 2:05pm
        
    end
    
    %Get the filetimes
    filetimes = [];
    filetimes = cellfun(@(x) getFileTime(x), filenames, 'UniformOutput', 0);
    filetimes = cell2mat(filetimes);
    
    %Based on which day you're processing, make index to the approapriate files
    ind = [];
    numMotifs = length(filt_audio);
    if i ~= 1
        %If it's not the day of lesion, just take the first and last 25 motifs
        ind(1,:) = 1:25;
        ind(2,:) = numMotifs-24:numMotifs;
        
        timeMarker = [i; i+.9]; %Day/night markers
    else
        %If it is the day of lesion, take 25 motif before and after lesion, then the last 25 of the day
        critInd = find(filetimes<critTime,1,'last');
        
        ind(1,:) = critInd-24:critInd;
        ind(2,:) = critInd+1:critInd+25;
        ind(3,:) = numMotifs-24:numMotifs;
        
        timeMarker = [i; i+.5; i+.9]; %Pre/post/night markers
    end
    
    %Do the calulation you care about
    for j = 1:size(ind,1)
        %Select the desired files
        subSet = filt_audio(ind(j,:));
        
        %Calculate total power in each motif spectrogram
        a = [];
%     a = cellfun(@(x) sum(sum(10.^x)), subSet, 'UniformOutput', 0);
        a = cellfun(@(x) sum(nanmean(windowTS(x.^2, 220, 44, 'pad', 'boxcar'),2)), subSet, 'UniformOutput', 0);
        aPower(j,:) = cell2mat(a);
    end
    
    %Shuffle data to the accumlating arrays
    runTime = [runTime; timeMarker];
    runPower = [runPower; aPower];
    aPower = []; timeMarker = [];
end

%Prep data for save
output.name = birdname;
output.time = runTime;
output.date = str2double(dates);
output.power = runPower;

m = exist([outputLocation, outName]);
if m ==2
    %File already exists
    clear('songPow');
    load([outputLocation, outName], 'songPow')
    songPow(end+1) = output;
else
    %No file yet created
    songPow = output;
end

%Save the updated data to file
save([outputLocation, outName], 'songPow')

display(['done with ' output.name])

% clear all





