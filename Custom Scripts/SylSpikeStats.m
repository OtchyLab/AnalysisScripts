function [Stats, CellStats]=SylSpikeStats(Location, Syls, StartRec, EndRec)
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

%If less than 4 argins, test all records in cell file
if nargin < 4
    StartRec = 0;
    EndRec = inf;
end

%If Syls is empty, set array to 1-10
if isempty(Syls)
    Syls = 1:10;
end

%Generate structures for all features to be analyzed
for i=1:length(Syls)
    Stats.(genvarname(['SpikeTimes' int2str(Syls(i))])) = [];
    Stats.(genvarname(['ISI' int2str(Syls(i))])) = [];
    Stats.(genvarname(['IFR' int2str(Syls(i))])) = [];
    Stats.(genvarname(['AveFR' int2str(Syls(i))])) = 0;
    Stats.(genvarname(['Bursts' int2str(Syls(i))])) = 0;
    Stats.(genvarname(['BurstFraction' int2str(Syls(i))])) = 0;
    Stats.(genvarname(['Sparseness' int2str(Syls(i))])) = 0;
    CellStats.(genvarname(['SpikeTimes' int2str(Syls(i))])) = [];
    CellStats.(genvarname(['SpikeTimesBreakout' int2str(Syls(i))])) = [];
    CellStats.(genvarname(['ISI' int2str(Syls(i))])) = [];
    CellStats.(genvarname(['IFR' int2str(Syls(i))])) = [];
    CellStats.(genvarname(['AveFR' int2str(Syls(i))])) = 0;
    CellStats.(genvarname(['Bursts' int2str(Syls(i))])) = 0;
    CellStats.(genvarname(['BurstFraction' int2str(Syls(i))])) = 0;
    CellStats.(genvarname(['Sparseness' int2str(Syls(i))])) = 0;
    CellStats.(genvarname(['SparsenessV' int2str(Syls(i))])) = [];
end
    

numfiles = length(presentCell);
%n = 1:numfiles;
SylLatency = 0.015; %***Set the amount of time (in sec) to collect spikes before Syl onset
q=1;

for i = 1:numfiles
    if (presentCell(i).filenum>=StartRec && presentCell(i).filenum<=EndRec) %&& (presentCell(i).filenum==elements{i}.filenum)
        for t = 1:length(elements)
            if elements{t}.filenum == presentCell(i).filenum
                EleInd = t;
            end
        end
        Stats(q).filename = presentCell(i).filenum;
        
        for SylInd = 1:length(Syls)
            %Find the times during which the bird is singing the selected
            %syllable; eliminate all other start/end times
            SylStarts = elements{EleInd}.segFileStartTimes(elements{EleInd}.segType(:)==Syls(SylInd));
            SylEnds = elements{EleInd}.segFileEndTimes(elements{EleInd}.segType(:)==Syls(SylInd));
        
            for j=1:length(SylStarts)
                %Calculate singing time (s) for the current syllable
                SylTime = SylEnds(j)-(SylStarts(j)-SylLatency); 
                Stats(q).(genvarname(['SingTime' int2str(Syls(SylInd))]))(j) = SylTime; %Calc total singing time for file
                
                %Get spikes for this syllable rep
                RepSpikes = presentCell(i).spikes(presentCell(i).spikes>=(SylStarts(j)-SylLatency) & presentCell(i).spikes<=SylEnds(j));
      
                %Check to see if there are any spikes for this section; if
                %not, fill with placeholder data

                %Store SpikeTimes for the Syl Rep, each aligned to sound onset
                if isempty(RepSpikes)
                    Stats(q).(genvarname(['SpikeTimes' int2str(Syls(SylInd))])) = [Stats(q).(genvarname(['SpikeTimes' int2str(Syls(SylInd))])); NaN];
                    CellStats.(genvarname(['SpikeTimes' int2str(Syls(SylInd))])) = [CellStats.(genvarname(['SpikeTimes' int2str(Syls(SylInd))])); NaN];
                else
                    Stats(q).(genvarname(['SpikeTimes' int2str(Syls(SylInd))])) = [Stats(q).(genvarname(['SpikeTimes' int2str(Syls(SylInd))])); RepSpikes-min(RepSpikes)];
                    CellStats.(genvarname(['SpikeTimes' int2str(Syls(SylInd))])) = [CellStats.(genvarname(['SpikeTimes' int2str(Syls(SylInd))])); RepSpikes-min(RepSpikes)];
                end
                
                %Calculate the ISI for spikes during the selected syllable
                k=1:length(RepSpikes)-1;
                RepISI = RepSpikes(k+1)-RepSpikes(k);
                Stats(q).(genvarname(['ISI' int2str(Syls(SylInd))])) = [Stats(q).(genvarname(['ISI' int2str(Syls(SylInd))])); RepISI];

                CellStats.(genvarname(['ISI' int2str(Syls(SylInd))])) = [CellStats.(genvarname(['ISI' int2str(Syls(SylInd))])); RepISI];

                %Calculate instantaneous firing rate during the selected syllable
                RepIFR = (RepISI).^-1;
                Stats(q).(genvarname(['IFR' int2str(Syls(SylInd))])) = [Stats(q).(genvarname(['IFR' int2str(Syls(SylInd))])); RepIFR];
                Stats(q).(genvarname(['AveFR' int2str(Syls(SylInd))])) = mean(Stats(q).(genvarname(['IFR' int2str(Syls(SylInd))])));
                
                
                CellStats.(genvarname(['IFR' int2str(Syls(SylInd))])) = [CellStats.(genvarname(['IFR' int2str(Syls(SylInd))])); RepIFR];

                %Calculate bursting spikes and bursting fraction
                k=1:length(RepIFR)-1;
                RepBurst = sum(RepIFR(k)>=150 & RepIFR(k+1)>=150);
                Stats(q).(genvarname(['Bursts' int2str(Syls(SylInd))])) = Stats(q).(genvarname(['Bursts' int2str(Syls(SylInd))])) + RepBurst;
                Stats(q).(genvarname(['BurstsFraction' int2str(Syls(SylInd))])) = Stats(q).(genvarname(['Bursts' int2str(Syls(SylInd))]))/length(Stats(q).(genvarname(['SpikeTimes' int2str(Syls(SylInd))])));
                
                CellStats.(genvarname(['Bursts' int2str(Syls(SylInd))])) = CellStats.(genvarname(['Bursts' int2str(Syls(SylInd))])) + RepBurst;
                
           end %End of reps for single syllable
        
        end %End of searching for a particular syllable
        q=q+1; %Index to next Stats record
        
    %elseif (presentCell(i).filenum~=elements{EleInd}.filenum)
    %    error(['Filenumber mismatch at i=' int2str(i)])
    end %End for checking the match of presentCell and elements
    
end %End of analyzing all files in presentCell

%Finish calculating cell-level features
for i = 1:length(Syls)
    CellStats.(genvarname(['AveFR' int2str(Syls(i))])) = mean(CellStats.(genvarname(['IFR' int2str(Syls(i))])));
    CellStats.(genvarname(['BurstFraction' int2str(Syls(i))])) = (CellStats.(genvarname(['Bursts' int2str(Syls(i))])))/length(CellStats.(genvarname(['IFR' int2str(Syls(i))])));
    
    index = find(CellStats.(genvarname(['SpikeTimes' int2str(Syls(i))])) == 0 | isnan(CellStats.(genvarname(['SpikeTimes' int2str(Syls(i))]))));

    %Separate out spike times for PSTH/Sparseness calcs
    for j = 1:length(index)
        if j~=length(index)
            CellStats.(genvarname(['SpikeTimesBreakout' int2str(Syls(i))])){j} = CellStats.(genvarname(['SpikeTimes' int2str(Syls(i))]))(index(j):max(index(j),index(j+1)-1));
        else
            CellStats.(genvarname(['SpikeTimesBreakout' int2str(Syls(i))])){j} = CellStats.(genvarname(['SpikeTimes' int2str(Syls(i))]))(index(j):max(index(j),end));
        end
    end
    Bins = 0:.003:max(CellStats.(genvarname(['SpikeTimes' int2str(Syls(i))])));%Use 3ms bins for PSTH
    psth = zeros(1,length(Bins));
    for k = 1:length(CellStats.(genvarname(['SpikeTimesBreakout' int2str(Syls(i))])))
        if (CellStats.(genvarname(['SpikeTimesBreakout' int2str(Syls(i))])){k} ==0 | isnan(CellStats.(genvarname(['SpikeTimesBreakout' int2str(Syls(i))])){k}))
            psth = psth + histc(CellStats.(genvarname(['SpikeTimesBreakout' int2str(Syls(i))])){k}, Bins);
        else
            psth = psth + histc(CellStats.(genvarname(['SpikeTimesBreakout' int2str(Syls(i))])){k}, Bins)';
        end
    end
    psth = psth/sum(psth);
    log_psth = zeros(1,length(psth));
    for m = 1:length(psth)
        if psth(m) == 0
            log_psth(m) = 0;
        else
            log_psth(m) = log(psth(m));
        end
    end
    CellStats.(genvarname(['Sparseness' int2str(Syls(i))])) = 1+sum(psth.*log_psth)/log(length(Bins));
    %CellStats.(genvarname(['SparsenessV' int2str(Syls(i))])) = cat(2,CellStats.(genvarname(['SparsenessV' int2str(Syls(i))])),CellStats.(genvarname(['Sparseness' int2str(Syls(i))])));
end
    


end %function end























