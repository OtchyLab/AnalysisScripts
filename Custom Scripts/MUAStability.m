%Generate summary figure for correlation/distance over days

%Constants (for this run...)
% dataName = 'C:\Users\Tim\Desktop\Temp Nif Data\Pur683\Pur683 HVC Activity (Ch4).mat';
% dataName = 'C:\Users\Tim\Desktop\Temp Nif Data\Pur690\Pur690 HVC Activity (Ch4).mat';
% dataName = 'C:\Users\Tim\Desktop\Temp Nif Data\Pur696\Pur696 HVC Activity (Ch3).mat';
% dataName = 'C:\Users\Tim\Desktop\Temp Nif Data\Pur755\Pur755 HVC Activity (Ch3).mat';

 dataName = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn046\Grn046 HVC Activity (Ch4).mat';
% dataName = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn121\Grn121 HVC Activity (Ch4).mat';
% dataName = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn141\Grn141 HVC Activity (Ch4).mat';

outputName =  'C:\Users\Tim\Desktop\Temp Nif Data\RecoveryData 1-6-2015.mat';
outputDir = 'C:\Users\Tim\Desktop\Temp Nif Data\';

win.start = 10;
win.end = 1;
lesion = true;

%If output doesn't exist, create it
% if ~exist('dataOut')
%     dataOut = [];
% end

%Load data from file
load(dataName, 'data');

%Grab the bird name and channel for later use
parts = regexp(dataName, '\', 'split');
birdName = parts{6};
birdChan = parts{end}(end-7:end-5);
figureName = parts{end}(1:end-4);


%Generate break indices
 [starts, ends, daybreaks, days] = genRendInd(data.files, win.start, win.end);
 %days(2) = 4;
 
 %Adjust baseline indices as needed
 if lesion
    baselineIdx = starts(1):data.lesionNum;
 else
    baselineIdx = starts(1):ends(1);
 end

%Calculate the recovery of mean HVC power
recoveryPow.m(1) = mean(mean(data.neuropowAbs(baselineIdx,:), 2));
recoveryPow.std(1) =std(mean(data.neuropowAbs(baselineIdx,:), 2));
for i = 2:length(days)
     recoveryPow.m(i) = mean(mean(data.neuropowAbs(starts(i):ends(i),:), 2));
     recoveryPow.std(i) = std(mean(data.neuropowAbs(starts(i):ends(i),:), 2));
end
recoveryPowNorm.m = recoveryPow.m./recoveryPow.m(1);
recoveryPowNorm.std = recoveryPow.std./recoveryPow.m(1);

h(1) = figure(1);
errorbar(days, recoveryPowNorm.m, recoveryPowNorm.std, '-xb')
title('Power Recovery')

%Calculate the recovery of correlation
[recoveryCorr, recoveryCorrNorm] = getBands(data.neuroCov, baselineIdx, starts, ends);

h(2) = figure(2);
errorbar(days, recoveryCorrNorm.m, recoveryCorrNorm.std, '-xb')
title('Correlation Recovery')

%Calculate the recovery of windowed correlation
[recoveryCorrWin, recoveryCorrWinNorm] = getBands(data.neuroCovWin, baselineIdx, starts, ends);

h(3) = figure(3);
errorbar(days, recoveryCorrWinNorm.m, recoveryCorrWinNorm.std, '-xb')
title('Windowed Correlation Recovery')

%Calculate the recovery of euclidean distance
[recoveryDist, recoveryDistNorm] = getBands(data.euclidDist, baselineIdx, starts, ends);

h(4) = figure(4);
errorbar(days, recoveryDistNorm.m, recoveryDistNorm.std, '-xb')
title('Recovered Distance')

%Calculate the recovery of windowed euclidean distance
[recoveryDistWin, recoveryDistWinNorm] = getBands(data.euclidDist_win, baselineIdx, starts, ends);

h(5) = figure(5);
errorbar(days, recoveryDistWinNorm.m, recoveryDistWinNorm.std, '-xb')
title('Recovered Windowed Distance')

%Send output to saved space
dataOut = [];
if exist(outputName, 'file')
    load(outputName,'dataOut')

    pntr = length(dataOut)+1;
else
    pntr = 1;
end
 
%ID stuff
dataOut(pntr).birdName = birdName;
dataOut(pntr).birdChan = birdChan;

%When and where
dataOut(pntr).days = days;
dataOut(pntr).win = win;

%Crunched data
dataOut(pntr).recoveryPow = recoveryPow;
dataOut(pntr).recoveryPowNorm = recoveryPowNorm;

dataOut(pntr).recoveryCorr = recoveryCorr;
dataOut(pntr).recoveryCorrNorm = recoveryCorrNorm;

dataOut(pntr).recoveryCorrWin = recoveryCorrWin;
dataOut(pntr). recoveryCorrWinNorm =  recoveryCorrWinNorm;

dataOut(pntr).recoveryDist = recoveryDist;
dataOut(pntr).recoveryDistNorm = recoveryDistNorm;

dataOut(pntr).recoveryDistWin = recoveryDistWin;
dataOut(pntr).recoveryDistWinNorm = recoveryDistWinNorm;

%Resave updated output data to file
save(outputName, 'dataOut');

%Save figures
savefig(h, [outputDir, figureName, '.fig'])

clear all

