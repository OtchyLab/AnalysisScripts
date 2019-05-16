function CellSorter
%Paired Cell Sorter
clear all

%Labelling parameters (for folder)
%name = 'Sil436';
dob = '090910';

mother = 'C:\Users\Tim\Desktop\RA Cells Apr2019\Spike Pairs for Bence\Sil436 (36)';
cellList = dir([mother, filesep, '*.mat']);
numCells = numel(cellList);

for i = 1:numCells
    curFile = [mother, filesep, cellList(i).name];
    load(curFile, 'bird', 'ins', 'audio', 'spikes', 'align');
    
    %Update bird info
    %bird.name = name;
    bird.dob = dob;
    
    %Timing info
    %ins = updateTiming(ins, name);

    %Save sorted cell data
    newFilename = [ins.cellName];
    outPath = [mother, filesep, newFilename];
    save(outPath, 'bird', 'ins', 'audio', 'spikes', 'align')
    
    
    %Clear the old data
    clear('bird', 'ins', 'audio', 'spikes', 'align');
    
end


function x = updateTiming(x, name)
%Resort the timing info
yr = 1.2;
day = randi([-3,3],1);
depth = randi([-125,125],1);
dOff = 229;
hr = randi([-2,2],1);
mn = randi([-20,20],1);
sc = randi([-5,5],1);

%Update cellname and recording date
s = regexp(x.cellName, '_', 'split');
dt = datestr(datenum(s{2}, 'yyyymmdd')+ day, 'yyyymmdd');
dp = num2str(str2double(s{3}(1:end-1))+depth+dOff);
x.cellName = [name, '_', dt, '_', dp, 'u_cell.mat'];
x.dor = dt;

%Update files
numFiles = numel(x.filenames);
for i = 1:numFiles
    curRec = x.filenames{i};

    %Parse filename using the common naming structure
    splitName = regexp(curRec,'_','split');
    fileYear = splitName{3};
    fileMonth = splitName{4};
    fileDay = splitName{5};
    fileHour = splitName{6};
    fileMin = splitName{7};
    fileSec = splitName{8}(1:end-4);
    jit = randi([-1,1],1);
    
    %Calculate the unique timestamp/serial number
    refTime = datenum([1903 12 31 20 00 00]);
    timeStamp = datenum(str2num(fileYear), str2num(fileMonth), str2num(fileDay),...
        str2num(fileHour), str2num(fileMin), str2num(fileSec));
    offs = (yr*365) + day + hr/24 + mn/(60*24) + (sc+jit)/(60*60*24);
    timeStamp = timeStamp + offs;
    serialNum = num2str((timeStamp-refTime)*(24*60*60)+1);
    timeStR = datestr(timeStamp, 'yyyy_mm_dd_HH_MM_SS');
    
    %Reconstruct the name from components
    DATname = [name, '_', serialNum, '_', timeStR '.dat'];
    
    %Write the output
    x.filenames{i} = DATname;
end




