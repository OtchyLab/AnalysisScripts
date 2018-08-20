function [Vals, Aves] = PlotUDStats(Location, Syls)
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
listD = who('CStats*D');
listU = [];

% for k = 1:length(listD)
%     listU(k) = regexprep(char(listD(k)), 'D', 'U');
% end

DAges = [117, 125, 125, 131];

% for i=1:length(Syls)
%     ValsU.(genvarname(['ISI' int2str(Syls(i))])) = [];
%     ValsU.(genvarname(['AveFR' int2str(Syls(i))])) = [];
%     ValsU.(genvarname(['BurstFraction' int2str(Syls(i))])) = [];
%     ValsU.(genvarname(['Sparseness' int2str(Syls(i))])) = [];
%     ValsU.(genvarname(['SpikeTimesBreakout' int2str(Syls(i))])) = [];
% end

for j=1:length(listD)
     AvesU.(genvarname(['AveFR' int2str(DAges(j))])) = [];
     AvesU.(genvarname(['BurstFraction' int2str(DAges(j))])) = [];
     AvesU.(genvarname(['Sparseness' int2str(DAges(j))])) = [];
     AvesU.(genvarname(['ISI' int2str(DAges(j))])) = [];
     AvesD.(genvarname(['AveFR' int2str(DAges(j))])) = [];
     AvesD.(genvarname(['BurstFraction' int2str(DAges(j))])) = [];
     AvesD.(genvarname(['Sparseness' int2str(DAges(j))])) = [];
     AvesD.(genvarname(['ISI' int2str(DAges(j))])) = [];
end

for i=1:length(Syls)
    ValsU.(genvarname(['ISI' int2str(Syls(i))])) = [];
    ValsU.(genvarname(['AveFR' int2str(Syls(i))])) = [];
    ValsU.(genvarname(['BurstFraction' int2str(Syls(i))])) = [];
    ValsU.(genvarname(['Sparseness' int2str(Syls(i))])) = [];
    ValsU.(genvarname(['SpikeTimesBreakout' int2str(Syls(i))])) = [];
    ValsD.(genvarname(['ISI' int2str(Syls(i))])) = [];
    ValsD.(genvarname(['AveFR' int2str(Syls(i))])) = [];
    ValsD.(genvarname(['BurstFraction' int2str(Syls(i))])) = [];
    ValsD.(genvarname(['Sparseness' int2str(Syls(i))])) = [];
    ValsD.(genvarname(['SpikeTimesBreakout' int2str(Syls(i))])) = [];
end

 for j=1:length(listD)
     CurFileD = char(listD(j));
     CurFileU = regexprep(char(listD(j)), 'D', 'U');
     for i=1:length(Syls)
        ValsD.(genvarname(['AveFR' int2str(Syls(i))])) = [ValsD.(genvarname(['AveFR' int2str(Syls(i))])); eval([CurFileD '.AveFR' int2str(Syls(i))])];
        AvesD.(genvarname(['AveFR' int2str(DAges(j))])) = [AvesD.(genvarname(['AveFR' int2str(DAges(j))])); eval([CurFileD '.AveFR' int2str(Syls(i))])];
        ValsD.(genvarname(['BurstFraction' int2str(Syls(i))])) = [ValsD.(genvarname(['BurstFraction' int2str(Syls(i))])); eval([CurFileD '.BurstFraction' int2str(Syls(i))])];
        AvesD.(genvarname(['BurstFraction' int2str(DAges(j))])) = [AvesD.(genvarname(['BurstFraction' int2str(DAges(j))])); eval([CurFileD '.BurstFraction' int2str(Syls(i))])];
        ValsD.(genvarname(['Sparseness' int2str(Syls(i))])) = [ValsD.(genvarname(['Sparseness' int2str(Syls(i))])); eval([CurFileD '.Sparseness' int2str(Syls(i))])];
        AvesD.(genvarname(['Sparseness' int2str(DAges(j))])) = [AvesD.(genvarname(['Sparseness' int2str(DAges(j))])); eval([CurFileD '.Sparseness' int2str(Syls(i))])];
        ValsD.(genvarname(['ISI' int2str(Syls(i))])) = [ValsD.(genvarname(['ISI' int2str(Syls(i))])); eval([CurFileD '.ISI' int2str(Syls(i))])];
        AvesD.(genvarname(['ISI' int2str(DAges(j))])) = [AvesD.(genvarname(['ISI' int2str(DAges(j))])); eval([CurFileD '.ISI' int2str(Syls(i))])];

        ValsU.(genvarname(['AveFR' int2str(Syls(i))])) = [ValsU.(genvarname(['AveFR' int2str(Syls(i))])); eval([CurFileU '.AveFR' int2str(Syls(i))])];
        AvesU.(genvarname(['AveFR' int2str(DAges(j))])) = [AvesU.(genvarname(['AveFR' int2str(DAges(j))])); eval([CurFileU '.AveFR' int2str(Syls(i))])];
        ValsU.(genvarname(['BurstFraction' int2str(Syls(i))])) = [ValsU.(genvarname(['BurstFraction' int2str(Syls(i))])); eval([CurFileU '.BurstFraction' int2str(Syls(i))])];
        AvesU.(genvarname(['BurstFraction' int2str(DAges(j))])) = [AvesU.(genvarname(['BurstFraction' int2str(DAges(j))])); eval([CurFileU '.BurstFraction' int2str(Syls(i))])];
        ValsU.(genvarname(['Sparseness' int2str(Syls(i))])) = [ValsU.(genvarname(['Sparseness' int2str(Syls(i))])); eval([CurFileU '.Sparseness' int2str(Syls(i))])];
        AvesU.(genvarname(['Sparseness' int2str(DAges(j))])) = [AvesU.(genvarname(['Sparseness' int2str(DAges(j))])); eval([CurFileU '.Sparseness' int2str(Syls(i))])];
        ValsU.(genvarname(['ISI' int2str(Syls(i))])) = [ValsU.(genvarname(['ISI' int2str(Syls(i))])); eval([CurFileU '.ISI' int2str(Syls(i))])];
        AvesU.(genvarname(['ISI' int2str(DAges(j))])) = [AvesU.(genvarname(['ISI' int2str(DAges(j))])); eval([CurFileU '.ISI' int2str(Syls(i))])];
     end
 end

 for j=1:length(unique(DAges))
     MeanAveFRD(j) = mean(AvesD.(genvarname(['AveFR' int2str(DAges(j))]))(~isnan(AvesD.(genvarname(['AveFR' int2str(DAges(j))])))));
     MeanBurstFractionD(j) = mean(AvesD.(genvarname(['BurstFraction' int2str(DAges(j))]))(~isnan(AvesD.(genvarname(['BurstFraction' int2str(DAges(j))])))));
     MeanSparsenessD(j) = mean(AvesD.(genvarname(['Sparseness' int2str(DAges(j))]))(~isnan(AvesD.(genvarname(['Sparseness' int2str(DAges(j))])))));
     
     MeanAveFRU(j) = mean(AvesU.(genvarname(['AveFR' int2str(DAges(j))]))(~isnan(AvesU.(genvarname(['AveFR' int2str(DAges(j))])))));
     MeanBurstFractionU(j) = mean(AvesU.(genvarname(['BurstFraction' int2str(DAges(j))]))(~isnan(AvesU.(genvarname(['BurstFraction' int2str(DAges(j))])))));
     MeanSparsenessU(j) = mean(AvesU.(genvarname(['Sparseness' int2str(DAges(j))]))(~isnan(AvesU.(genvarname(['Sparseness' int2str(DAges(j))])))));
 end
     
Colors = ['b'; 'm'; 'c'; 'r'; 'g'; 'b'; 'r'; 'k'];
figure
for i=1:length(Syls)
    subplot(2,4,i)
    scatter(ValsU.(genvarname(['AveFR' int2str(Syls(i))])), ValsD.(genvarname(['AveFR' int2str(Syls(i))])), Colors(i));
    xlabel('Average Firing Rate U (Hz)');
    ylabel('Average Firing Rate D (Hz)');
    title(['Average Firing Rate for Syl ' int2str(Syls(i)) ' (U vs D)']);
    xlim([0, 500]);
    ylim([0, 500]);
    hold on
    plot([0, 500], [0, 500], 'k')
    hold on
end

figure
for i=1:length(Syls)
    subplot(2,4,i)
    scatter(ValsU.(genvarname(['Sparseness' int2str(Syls(i))])), ValsD.(genvarname(['Sparseness' int2str(Syls(i))])), Colors(i));
    xlabel('Sparseness U');
    ylabel('Sparseness D');
    title(['Sparseness for Syl ' int2str(Syls(i)) ' (U vs D)']);
    xlim([0, 1]);
    ylim([0, 1]);
    hold on
    plot([0, 1], [0, 1], 'k')
    hold on
end

figure
for i=1:length(Syls)
    subplot(2,4,i)
    scatter(ValsU.(genvarname(['BurstFraction' int2str(Syls(i))])), ValsD.(genvarname(['BurstFraction' int2str(Syls(i))])), Colors(i));
    xlabel('Bursting Fraction U');
    ylabel('Bursting Fraction D');
    title(['Burst Fraction for Syl ' int2str(Syls(i)) ' (U vs D)']);
    xlim([0, 1]);
    ylim([0, 1]);
    hold on
    plot([0, 1], [0, 1], 'k')
    hold on
end

figure
subplot(3,1,1)
scatter(MeanAveFRU, MeanAveFRD, Colors(1));
xlabel('Undirected Average Firing Rate (Hz)');
ylabel('DirectedAverage Firing Rate (Hz)');
title(['Average Firing Rate (Undirected vs Directed)']);
xlim([0, 500]);
ylim([0, 500]);
hold on
plot([0, 500], [0, 500], 'k')
hold on

subplot(3,1,2)
scatter(MeanBurstFractionU, MeanBurstFractionD, Colors(2));
xlabel('Undirected Bursting Fraction');
ylabel('Directed Bursting Fraction');
title(['Bursting Fraction (Undirected vs Directed)']);
xlim([0, 1]);
ylim([0, 1]);
hold on
plot([0, 1], [0, 1], 'k')
hold on

subplot(3,1,3)
scatter(MeanSparsenessU, MeanSparsenessD, Colors(3));
xlabel('Undirected Sparseness');
ylabel('Directed Sparseness');
title(['Sparseness (Undirected vs Directed)']);
xlim([0, 1]);
ylim([0, 1]);
hold on
plot([0, 1], [0, 1], 'k')
hold on


a=1







end %Function End














