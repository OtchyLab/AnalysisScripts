function [Vals, Aves] = PlotStats(Location, Syls)
%Function calculates both File and Cell statistics for the provided .mat
%files.  Both the 'cell' and 'annotation' files are needed for this
%function.
%NOTE: this function looks at all spikes produced during times that a 
%specified syllable occurs.  There is no distinction between motifs.
% Input: 
% Location - absolute path of the times, cell, and annotation files
%           (string)
% Syls - an array specifying which syllables to calc statistics for; an
%           empty set ([]) is a short cut for all numbered 1-10
% StartRec - the first record number (given in the .filenum field)
%           to be processed (remove to process all)
% EndRec - the last record number (given in the .filenum field)
%           to be processed (remove to process all)
%--------------------------------------------------------------------------
% Output:
% Stats - Spike stats sorted by file
% CellStats - Spike stats for the entire set of selected recordings
%--------------------------------------------------------------------------

%Comment out below to use command line Location
Location = 'C:\Users\Tim\Desktop\SpikeStats Test Folder';
%Location = '/Users/MacDaddy/Desktop/SpikeStats Test Folder';

%Check to see if Location exists; if so, set as CD, index, and load files
if ~isdir(Location)
    error('Location does not exist')
    return
end
cd(Location);

%Load and verify *cell.mat
files = dir('Stats*.mat');
load(files(1).name); %***Only looks at the first cell-file found***
list = who('CStats*U');

DAges = [111, 111, 112, 112, 112, 114, 117, 120, 125, 125, 131];

for i=1:length(Syls)
    Vals.(genvarname(['ISI' int2str(Syls(i))])) = [];
    Vals.(genvarname(['AveFR' int2str(Syls(i))])) = [];
    Vals.(genvarname(['BurstFraction' int2str(Syls(i))])) = [];
    Vals.(genvarname(['Sparseness' int2str(Syls(i))])) = [];
    Vals.(genvarname(['SpikeTimesBreakout' int2str(Syls(i))])) = [];
end

for j=1:length(list)
     Aves.(genvarname(['AveFR' int2str(DAges(j))])) = [];
     Aves.(genvarname(['BurstFraction' int2str(DAges(j))])) = [];
     Aves.(genvarname(['Sparseness' int2str(DAges(j))])) = [];
     Aves.(genvarname(['ISI' int2str(DAges(j))])) = [];
end

 for j=1:length(list)
     CurFile = char(list(j));
     for i=1:length(Syls)
        Vals.(genvarname(['AveFR' int2str(Syls(i))])) = [Vals.(genvarname(['AveFR' int2str(Syls(i))])); eval([CurFile '.AveFR' int2str(Syls(i))])];
        Aves.(genvarname(['AveFR' int2str(DAges(j))])) = [Aves.(genvarname(['AveFR' int2str(DAges(j))])); eval([CurFile '.AveFR' int2str(Syls(i))])];
        Vals.(genvarname(['BurstFraction' int2str(Syls(i))])) = [Vals.(genvarname(['BurstFraction' int2str(Syls(i))])); eval([CurFile '.BurstFraction' int2str(Syls(i))])];
        Aves.(genvarname(['BurstFraction' int2str(DAges(j))])) = [Aves.(genvarname(['BurstFraction' int2str(DAges(j))])); eval([CurFile '.BurstFraction' int2str(Syls(i))])];
        Vals.(genvarname(['Sparseness' int2str(Syls(i))])) = [Vals.(genvarname(['Sparseness' int2str(Syls(i))])); eval([CurFile '.Sparseness' int2str(Syls(i))])];
        Aves.(genvarname(['Sparseness' int2str(DAges(j))])) = [Aves.(genvarname(['Sparseness' int2str(DAges(j))])); eval([CurFile '.Sparseness' int2str(Syls(i))])];
        Vals.(genvarname(['ISI' int2str(Syls(i))])) = [Vals.(genvarname(['ISI' int2str(Syls(i))])); eval([CurFile '.ISI' int2str(Syls(i))])];
        Aves.(genvarname(['ISI' int2str(DAges(j))])) = [Aves.(genvarname(['ISI' int2str(DAges(j))])); eval([CurFile '.ISI' int2str(Syls(i))])];
        %        Vals.(genvarname(['SpikeTimesBreakout' int2str(Syls(i))])) =
%        [Vals.(genvarname(['SpikeTimesBreakout' int2str(Syls(i))])); eval([CurFile '.SpikeTimesBreakout' int2str(Syls(i))])];
     end
%      MeanAveFR(j) = mean(Aves.(genvarname(['AveFR' int2str(DAges(j))]))(~isnan(Aves.(genvarname(['AveFR' int2str(DAges(j))])))));
%      MeanBurstFraction(j) = mean(Aves.(genvarname(['BurstFraction' int2str(DAges(j))]))(~isnan(Aves.(genvarname(['BurstFraction' int2str(DAges(j))])))));
%      MeanSparseness(j) = mean(Aves.(genvarname(['Sparseness' int2str(DAges(j))]))(~isnan(Aves.(genvarname(['Sparseness' int2str(DAges(j))])))));
 end

 for j=1:length(unique(DAges))
     MeanAveFR(j) = mean(Aves.(genvarname(['AveFR' int2str(DAges(j))]))(~isnan(Aves.(genvarname(['AveFR' int2str(DAges(j))])))));
     MeanBurstFraction(j) = mean(Aves.(genvarname(['BurstFraction' int2str(DAges(j))]))(~isnan(Aves.(genvarname(['BurstFraction' int2str(DAges(j))])))));
     MeanSparseness(j) = mean(Aves.(genvarname(['Sparseness' int2str(DAges(j))]))(~isnan(Aves.(genvarname(['Sparseness' int2str(DAges(j))])))));
 end
     
Colors = ['y'; 'm'; 'c'; 'r'; 'g'; 'b'; 'y'; 'k'];
figure
for i=1:length(Syls)
    subplot(2,4,i)
    scatter(DAges, Vals.(genvarname(['AveFR' int2str(Syls(i))])), Colors(i));
    xlabel('Age (dph)');
    ylabel('Average Firing Rate (Hz)');
    title(['Average Firing Rate for Syl ' int2str(Syls(i)) ' vs Age']);
    xlim([105, 140]);
    ylim([0, 500]);
    hold on
end

figure
for i=1:length(Syls)
    subplot(2,4,i)
    scatter(DAges, Vals.(genvarname(['Sparseness' int2str(Syls(i))])), Colors(i));
    xlabel('Age (dph)');
    ylabel('Sparseness');
    title(['Sparseness for Syl ' int2str(Syls(i)) ' vs Age']);
    xlim([105, 140]);
    ylim([0, 1]);
    hold on
end

figure
for i=1:length(Syls)
    subplot(2,4,i)
    scatter(DAges, Vals.(genvarname(['BurstFraction' int2str(Syls(i))])), Colors(i));
    xlabel('Age (dph)');
    ylabel('Bursting Fraction');
    title(['Burst Fraction for Syl ' int2str(Syls(i)) ' vs Age']);
    xlim([105, 140]);
    ylim([0, 1]);
    hold on
end

figure
scatter(unique(DAges), MeanAveFR, Colors(1));
xlabel('Age (dph)');
ylabel('Average Firing Rate (Hz)');
title(['Average Firing Rate vs Age']);
xlim([105, 140]);
ylim([0, 500]);
hold on

fitline = polyfit(unique(DAges), MeanAveFR,1); % least squares fitting to a line
a1 = fitline(2); % y-intercept of the fitted line
a2 = fitline(1); % slope of fitted lines
fit = a2*unique(DAges)+a1;
plot(unique(DAges),fit,'k')

figure
scatter(unique(DAges), MeanBurstFraction, Colors(2));
xlabel('Age (dph)');
ylabel('Bursting Fraction');
title(['Burst Fraction for Syl vs Age']);
xlim([105, 140]);
ylim([0, 1]);
hold on

fitline = polyfit(unique(DAges), MeanBurstFraction,1); % least squares fitting to a line
a1 = fitline(2); % y-intercept of the fitted line
a2 = fitline(1); % slope of fitted lines
fit = a2*unique(DAges)+a1;
plot(unique(DAges),fit,'k')

figure
scatter(unique(DAges), MeanSparseness, Colors(3));
xlabel('Age (dph)');
ylabel('Sparseness');
title(['Sparseness for Syl vs Age']);
xlim([105, 140]);
ylim([0, 1]);
hold on

fitline = polyfit(unique(DAges), MeanSparseness,1); % least squares fitting to a line
a1 = fitline(2); % y-intercept of the fitted line
a2 = fitline(1); % slope of fitted lines
fit = a2*unique(DAges)+a1;
plot(unique(DAges),fit,'k')




a=1







end %Function End














