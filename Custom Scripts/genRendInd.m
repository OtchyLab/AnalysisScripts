function [starts, ends, daybreaks, days] = genRendInd(filelist, startTime, endTime)
% 
%
%
%

%Cycle through each file and extract the date/time info
dateNums = [];
timeNums = [];
for i = 1:length(filelist)
    dateNums(i) = getFileTime(filelist{i});
    
    strparts = regexp(filelist{i},'_', 'split');
    th = str2double(strparts{6});
    tm = str2double(strparts{7});
    timeNums(i) = th + (tm/60);
end

days = (floor(dateNums-floor(dateNums(1))))+1;
daybreaks = [1,  find(diff(days) ~= 0) + 1, length(filelist)]; %Index number of first rendition of the day
days = unique(days); %day number that each row corresponds to

for i = 1:(length(daybreaks)-1)
    above = find(timeNums(daybreaks(i):(daybreaks(i+1)-1)) >= startTime, 1, 'first');
    if ~isempty(above)
        starts(i) = daybreaks(i) + above - 1;
    else
        starts(i) = NaN;
    end
    
    below = find(timeNums(daybreaks(i):(daybreaks(i+1)-1)) <= endTime, 1, 'last');
    if ~isempty(below)
        ends(i) = daybreaks(i) + below - 1;
    else
        ends(i) = NaN;
    end
    
    if starts(i) > ends(i)
        starts(i) = NaN;
        ends(i) = NaN;
    end
end

daybreaks = daybreaks(1:(end-1)); 