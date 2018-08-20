function [spkTrigAmp, spkTrigAM, spkTrigFM, spkTrigEnt, spkTrigPeakFreq, spkTrigPitch, spkTrigPitchGoodness, numSpikes] = prelimSpikeTriggeredFeatureAverages(cell , winSec, motifs)

%NOTE: no warping, features should just be computed at 40000Hz, features 
%should be checked for matching with SAP and pitch should be computed in 
%alternate fashions. 

if(~exist('motifs'))
    motifs = [1:length(cell.motifData)];
end

sampPerWin = 409;
advWin = 60; %number of 44.1 samples the feature estimator is advanced at each time.
noiseRatio = .5;
[tapers, concs] = dpss(sampPerWin, 1.5);

win = round(winSec*44100);
win = round(win./advWin);

%The average
numSpikes = 0;
nNumSkipped = 0;
spkTrigAmp = zeros(win*2+1,1);
spkTrigAM = zeros(win*2+1,1);
spkTrigFM = zeros(win*2+1,1);
spkTrigEnt = zeros(win*2+1,1);
spkTrigPeakFreq = zeros(win*2+1,1);
spkTrigPitch = zeros(win*2+1,1);
spkTrigPitchGoodness = zeros(win*2+1,1);

varAmp = 0;
varAM = 0;
varFM = 0;
varEnt = 0;
varPeakFreq = 0;
varPitch = 0;
varPitchGoodness = 0;

numMotifs = 0;
for(nMotif=motifs)
    motif = cell.motifData{nMotif};
    if(isfield(motif, 'spikeNdx'))
        audio = motif.audio;

        %This function automatically resamples to 44.1 (for now atleast to
        %match SAP.
        [specDeriv, fAmp, fAM, fFM, fEntropy, fPeakFreq, fPitch, fPitchGoodness] = SAP_computeLocalFeatures(audio, cell.exper.desiredInSampRate, sampPerWin, advWin, noiseRatio, tapers);
        
        numMotifs = numMotifs + 1;
        varAmp = varAmp + var(fAmp);
        varAM = varAmp + var(fAM);
        varFM = varAmp + var(fFM);
        varEnt = varAmp + var(fEntropy);
        varPeakFreq = varAmp + var(fPeakFreq);
        varPitch = varAmp + var(fPitch);
        varPitchGoodness = varAmp + var(fPitchGoodness);
        
        %Spike Ndx need to be scaled to be at 44.1, and then mapped to the
        %indices of spectral derivative space.
        spikeNdx = round(motif.spikeNdx .* (44100/40000));
        spikeNdx = round(spikeNdx ./ advWin);

        %Compute spike triggered average the slow way?
        for(nSpk = 1:length(spikeNdx))       
            spkNdx = spikeNdx(nSpk);
            startWin = spkNdx - win;
            endWin = spkNdx + win;
            if(startWin>=1 & endWin<=length(fAmp))
                numSpikes = numSpikes + 1;
                spkTrigAmp = spkTrigAmp + fAmp(startWin:endWin)'; amp(numSpikes,:) = fAmp(startWin:endWin)';
                spkTrigAM = spkTrigAM + fAM(startWin:endWin)'; am(numSpikes,:) = fAM(startWin:endWin)';
                spkTrigFM = spkTrigFM + fFM(startWin:endWin)'; fm(numSpikes,:) = fFM(startWin:endWin)';
                spkTrigEnt = spkTrigEnt + fEntropy(startWin:endWin)'; ent(numSpikes,:) = fEntropy(startWin:endWin)';
                spkTrigPeakFreq = spkTrigPeakFreq + fPeakFreq(startWin:endWin)'; peakfreq(numSpikes,:) = fPeakFreq(startWin:endWin)';
                spkTrigPitch = spkTrigPitch + fPitch(startWin:endWin)'; pitch(numSpikes,:) = fPitch(startWin:endWin)';
                spkTrigPitchGoodness = spkTrigPitchGoodness + fPitchGoodness(startWin:endWin)'; pitchgood(numSpikes,:) = fPitchGoodness(startWin:endWin)';  
            else
                nNumSkipped = nNumSkipped + 1;
            end
        end   
    end
end

if(numSpikes > 0)
    spkTrigAmp = spkTrigAmp / numSpikes;
    spkTrigAM = spkTrigAM / numSpikes;
    spkTrigFM = spkTrigFM / numSpikes;
    spkTrigEnt = spkTrigEnt / numSpikes;
    spkTrigPeakFreq = spkTrigPeakFreq / numSpikes;
    spkTrigPitch = spkTrigPitch / numSpikes;
    spkTrigPitchGoodness = spkTrigPitchGoodness / numSpikes;    
else
    error('No spikes to average');
end

varAmp / numMotifs
varAmp / numMotifs
varAmp / numMotifs
varAmp / numMotifs
varAmp / numMotifs
varAmp / numMotifs
varAmp / numMotifs
  
if(nNumSkipped > 0)
    warning(['SpikeTrigAvg: skipped ', num2str(nNumSkipped), 'spikes due to truncation.']);
end
    