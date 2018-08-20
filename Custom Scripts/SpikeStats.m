function [Stats, CellStats]=SpikeStats(Location, StartRec, EndRec)
%Function calculates both File and Cell statistics for the provided .mat
%files.  Both the 'cell' and 'annotation' files are needed for this
%function.
%NOTE: this function looks at all spikes produced during times that singing 
%(and not calls or un-ID'ed syllables) occurs.  There is no distinction
%between syllable types or motifs.
% Input: 
% Location - absolute path of the times, cell, and annotation files
%           (string)
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
files = dir('*cell.mat');
load(files(1).name); %***Only looks at the first cell-file found***
if ~isfield(presentCell, 'spikes')
    error('presentCell structure not loaded/formatted correctly')
    return
end

 %Load and verify *annotation.mat
files = dir('*annotation.mat');
load(files(1).name); %***Only looks at the first annotation-file found***
if ~exist('elements')
    error('elements array not loaded/formatted correctly')
    return
end

%If less than 3 argins, test all records in cell file
if nargin < 3
    StartRec = 0;
    EndRec = inf;
end

numfiles = min(length(presentCell), length(elements));
n = 1:numfiles;
SongStarts = [];
SongEnds = [];
CellStats.ISI = []; CellStats.IFR = []; CellStats.Bursts = 0;
q = 1;

for i = 1:numfiles
    if (presentCell(i).filenum>=StartRec & presentCell(i).filenum<=EndRec) & (presentCell(i).filenum==elements{i}.filenum)
        Stats(q).filename = presentCell(i).filenum;
        %Find the times during which the bird is singing
        %Eliminate calls and unknown syllables
        segFileStartTimes = elements{i}.segFileStartTimes(elements{i}.segType(:)<100);
        segFileEndTimes = elements{i}.segFileEndTimes(elements{i}.segType(:)<100);

        %Initialize structures for stripping out song times
        CurSyl = 1; Bouts = 1;
        SongStarts(Bouts) = segFileStartTimes(1);
        SongEnds(Bouts) = segFileEndTimes(1);

        %If gaps between syllables are less that 100ms, the bout is poly-syllabic
        while CurSyl<length(segFileStartTimes)
            if (segFileStartTimes(CurSyl+1)-segFileEndTimes(CurSyl)) <= 0.1
                SongEnds(Bouts) = segFileEndTimes(CurSyl+1);
            else
                Bouts = Bouts +1;
                SongStarts(Bouts) = segFileStartTimes(CurSyl+1);
                SongEnds(Bouts) = segFileEndTimes(CurSyl+1);
            end    
            CurSyl=CurSyl+1;
        end

        SpikeTimes = []; SylLatency = 0.05; BoutBurst =[];
        Stats(q).ISI = []; Stats(q).IFR = []; Stats(q).Bursts = 0; Stats(q).BurstFraction = 0; Stats(q).SingTime = 0;
        
        for j=1:Bouts
            BoutSpikes = presentCell(i).spikes(presentCell(i).spikes>=(SongStarts(j)-SylLatency) & presentCell(i).spikes<=SongEnds(j));
            %Calculate singing time (s) for bout
            SingTime = SongEnds(j)-(SongStarts(j)-SylLatency); 
            Stats(q).SingTime = Stats(q).SingTime + SingTime; %Calc total singing time for file
            SpikeTimes = [SpikeTimes; BoutSpikes]; %Build file of all file spikes during song

            %Calculate the ISI for spikes in this file
            k=1:length(BoutSpikes)-1;
            BoutISI = BoutSpikes(k+1)-BoutSpikes(k);
            Stats(q).ISI = [Stats(q).ISI; BoutISI];
            CellStats.ISI = [CellStats.ISI; BoutISI];

            %Calculate instantaneous firing rate
            BoutIFR = (BoutISI).^-1;
            Stats(q).IFR = [Stats(q).IFR; BoutIFR];
            CellStats.IFR = [CellStats.IFR; BoutIFR];
            Stats(q).AveFR = mean(Stats(q).IFR);

            %Calculate bursting spikes and bursting fraction
            k=1:length(BoutIFR)-1;
            BoutBurst = sum(BoutIFR(k)>=150 & BoutIFR(k+1)>=150);
            Stats(q).Bursts = Stats(q).Bursts + BoutBurst;
            CellStats.Bursts = CellStats.Bursts + BoutBurst;
            Stats(q).BurstFraction = (Stats(q).Bursts/length(Stats(q).IFR));
        end
        %Calculate sparseness file
        SparseBins = max(min(SpikeTimes)-.1, 0):.003:max(SpikeTimes); %Use 3ms bins over the length of the file
        psth = 0;
        log_psth = 0;
        SparsenessV = [];
        %for m=1:length(SpikeTimes)
            psth = psth + histc(SpikeTimes,SparseBins)';
        %end
        psth = psth/sum(psth); %Normalize PSTH
        for m = 1:length(psth)
            if psth(m) == 0;
                log_psth(m) = 0;
            else
                log_psth(m) = log(psth(m));
            end
        end
        Stats(q).Sparseness = 1+sum(psth.*log_psth)/log(length(SparseBins));
        SparsenessV = [SparsenessV; Stats(q).Sparseness];
        q=q+1;
    elseif (presentCell(i).filenum~=elements{i}.filenum)
        error(['Filenumber mismatch at i=' int2str(i)])
    end
end
CellStats.AveIFR = mean(CellStats.IFR);
CellStats.BurstFraction = (CellStats.Bursts/length(CellStats.IFR));
end %function end