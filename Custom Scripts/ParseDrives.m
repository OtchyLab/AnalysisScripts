%This script is to pull the relevant data out of Farhan's "Drive" files
%with pitch and timing data.

%Pull the drive into the local variable
drive = drive4;

%Parse out the data into useful arrays
renditions = size(drive,1);
for i = 1:renditions
   date(i) = datenum(drive{i}(1:10)); % date
   hour(i) = str2num(drive{i}(12:13)); % hour
   minute(i) = str2num(drive{i}(15:16)); % minute
   %if str2num(drive{i}(26))<=1
   %    pitchEst(i) = str2num(drive{i}(26:31)); % pitch Estimate
   %else
    pitchEst(i) = str2num(drive{i}(26:30)); % pitch Estimate
   %end
end

%Find unique dates
uDates = sort(unique(date))';

%Calculate period means
for i = 1:size(uDates,1)
   pointDate{i} = datestr(uDates(i),'mm/dd/yy');
   pitchMean(i) = mean(pitchEst(date(:)==uDates(i) & hour(:)==11));
   pitchStd(i) = std(pitchEst(date(:)==uDates(i) & hour(:)==11));
end

pointDate = pointDate';
pitchMean = pitchMean';
pitchStd = pitchStd';