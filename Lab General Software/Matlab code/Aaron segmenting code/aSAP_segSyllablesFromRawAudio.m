function [syllStartTimes, syllEndTimes, noiseEst, noiseStd, soundEst, soundStd] = aSAP_segSyllablesFromRawAudio(audio, fs)
%presume audio is filtered.

bDebug = false;

syllStartTimes = [];
syllEndTimes = [];
filesForManual = [];

fMinSyllDuration = .007; %secs
fMinIntervalDuration = .007; %secs
fMaxSyllDuration = 1; %sec

%utilizes several thresholds and criteria
%1. syllalbe trig threshold power - the syllable has to reach this
%   threshold at least once.
%2. syllable cont threshold power - all time between this power is
%   considered the syllable.
%3. min interval duration.
%4. min syllable duration.

%5. max syllable duration

%demean the audio.
audio = audio - mean(audio);

%include only very relevant power:
%below 8600Hz
%filt2.order = 50; %sufficient for 44100Hz of lower
%filt2.win = hann(filt2.order+1);
%filt2.cutoff = 8600; %Hz
%filt2.fs = fs;
%filt2.lpf = fir1(filt2.order, filt2.cutoff/(filt2.fs/2), 'low', filt2.win);
%audio = filtfilt(filt2.lpf, 1, audio);

%above 860Hz
filt3.order = 50; %sufficient for 44100Hz of lower
filt3.win = hann(filt3.order+1);
filt3.cutoff = 860; %Hz
filt3.fs = fs;
filt3.hpf = fir1(filt3.order, filt3.cutoff/(filt3.fs/2), 'high', filt3.win);
audio = filtfilt(filt3.hpf, 1, audio);

%compute power
audioPow= audio.^2; 

%smooth the power, lpf:
filt.order = 100; 
filt.win = hann(filt.order+1);
filt.cutoff = 50; %Hz
filt.fs = fs;
filt.lpf = fir1(filt.order, filt.cutoff/(filt.fs/2), 'low', filt.win);
audioPow = filtfilt(filt.lpf, 1, audioPow);

%compute log Pow
audioLogPow = log(audioPow + eps);


%estimate the noise and sound levels by finding the lowest and highest peaks in the estimated probability
%distribution:
%[f,x] = ksdensity(audioLogPow); %estimate probability distribtuion.
%f = f.*(f>.01);
%df = diff(f);
%peaks = find(df(1:end-1)>0 & df(2:end) < 0); %we assume the first peak is the noise peak
%troughs = find(df(1:end-1)<0 & df(2:end) > 0);
%if(length(peaks) <= 1 || length(troughs) == 0)
%    return;
%end
%noiseEst = x(peaks(1));  
%classEst = x(troughs((find(x(troughs)>noiseEst))));  %the first trough after the noise represents a divider.
%if(length(classEst) == 0)
%    return;
%else
%    classEst = classEst(1);
%end
%soundEst = x(peaks(end));
%noiseVar = var(audioLogPow(audioLogPow<classEst));
%noiseMean = mean(audioLogPow(audioLogPow<classEst));
%Set sound thresholds based on these estimates:
%thresTrig = classEst + .3 * (soundEst - classEst);
%thresSyll = classEst + .3 * (noiseEst - classEst);

%Estimate threshold to discrimate between sound and noise using a mixture
%of two gaussians model.
[noiseEst, soundEst, noiseStd, soundStd] = aSAP_estimateTwoMeans(audioLogPow);    
if(noiseEst>soundEst)
    return;
end
%Compute the optimal classifier between the two gaussians...
p(1) = 1/(2*soundStd^2) - 1/(2*noiseStd^2);
p(2) = (noiseEst)/(noiseStd^2) - (soundEst)/(soundStd^2);
p(3) = (soundEst^2)/(2*soundStd^2) - (noiseEst^2)/(2*noiseStd^2) + log(soundStd/noiseStd);
disc = roots(p);
disc = disc(find(disc>noiseEst & disc<soundEst));
if(length(disc)==0)
    return;
end
disc = disc(1);

%Set the thresholds based on these estimates
thresSyll = disc + 0.2 * (soundEst-disc);
thresTrig = soundEst - 0.4 * (soundEst - disc);

if(bDebug)
    figure(1115);
    clf;
    s1 = subplot(2,1,1);
    plot([0:length(audio)-1]/fs,audioLogPow);
    line(xlim,[noiseEst,noiseEst], 'Color', 'red');
    line(xlim,[soundEst,soundEst], 'Color', 'green');
    line(xlim,[thresTrig,thresTrig], 'Color', 'black');
    line(xlim,[thresSyll,thresSyll], 'Color', 'blue');
end

%Find threshold crossings:
[trigCross, junk] = detectThresholdCrossings(audioLogPow, thresTrig, true);
[syllUpCross, syllDownCross] = detectThresholdCrossings(audioLogPow, thresSyll, true);

if(length(syllUpCross) > 0 | length(syllDownCross) > 0)
    %Eliminated extraneous end crossing...
    if(syllUpCross(1) == 1)
        syllUpCross = syllUpCross(2:end);
    end
    if(syllDownCross(end) == length(audioLogPow))
        syllDownCross = syllDownCross(1:end-1);
    end

    %Determine syllables present
    nSyll = 0;
    beginSyll = [];
    endSyll = [];
    for(nTrig = 1:length(trigCross))
        up = find(syllUpCross<trigCross(nTrig));
        down = find(syllDownCross>trigCross(nTrig));
        if((length(up)>0) & (length(down)>0))
            nSyll = nSyll + 1;
            beginSyll(nSyll) = syllUpCross(up(end));
            endSyll(nSyll) = syllDownCross(down(1));
        end
    end

    if(length(beginSyll) > 2)
        %Remove small intervals
        intervals = (beginSyll(2:end) - endSyll(1:end-1)) ./ fs;
        realGapNdx = find(intervals > fMinIntervalDuration);
        beginSyll = beginSyll([1,realGapNdx+1]);
        endSyll = endSyll([realGapNdx,length(endSyll)]);
    end
        
    %Remove syllables that are too short or long
    durations = (endSyll - beginSyll) / fs;
    realSyll = find((durations > fMinSyllDuration) & (durations < fMaxSyllDuration));
    beginSyll = beginSyll(realSyll);
    endSyll = endSyll(realSyll);

    syllStartTimes = (beginSyll -1)/ fs;
    syllEndTimes = (endSyll-1) / fs;
end

if(bDebug)
    s2 = subplot(2,1,2);
    [SAP_Feats, m_spec_deriv] = aSAP_generateASAPFeatures(audio, fs);
    aSAP_displaySpectralDerivative(m_spec_deriv, Parameters);

    allCross = [syllStartTimes, syllEndTimes];
    for(i = 1:length(allCross))
        axes(s1);
        line([allCross(i), allCross(i)], ylim, 'Color', 'red');
        axes(s2);
        line([allCross(i), allCross(i)], ylim, 'Color', 'red');
    end
    linkaxes([s1,s2], 'x')
    figure(gcf);
    pause;
end





