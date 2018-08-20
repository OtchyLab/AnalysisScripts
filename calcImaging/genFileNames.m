% Script to get naming data from csv files for .wav/.mov files
% navigate to the directory containing the .mov files
cd '/Users/sadiela/Desktop/LW9'
% pull out list of all csv files in the directory
mov_files = dir('*.mov');
% set target directory
targetDirectory = "/Users/sadiela/Desktop/calc_data";

numFiles = length(mov_files);
for i = 1:numFiles
    curFileName = mov_files(i).name;
    
    % load data from csv file... all necessary info is in the .mov files so
    % unnecessary???
    %data = csv2cell(curFileName, 'fromfile');
    %date_time = curFileName;
    
    % get birdname from folder in which the current file is located
    birdName = getBirdName(pwd);
    
    % YYYY-MM-DD
    
    % get date from string date/time
    [Y, M, D, H, MN, S] = getDateVec(curFileName);
    
    % create serial time stamp
    refTime = datenum([1903 12 31 20 00 00]);
    timeStamp = datenum([Y, M, D, H, MN, S]);
    serialNum = num2str((timeStamp-refTime)*(24*60*60)+1);
    
    % create file name
    nameFile(i).name = formFileName(birdName, serialNum, Y, M, D, H, MN, S);
end
% now I have a struct with all the future names of the .wav and .mov files



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 x = 1:numel(raw);
    v = raw';

    fs = rate_new; %desired frequency
    n = round((numel(raw)/rate_orig) * fs);

    xq = linspace(1, numel(v), n);

    audioData = interp1(x, v, xq);



















