function filteraudio(directory)
%fs is the sampling rate.

%syllStartTimes: the start time of each syllable.  This is not the onset of
%sound, but the time at which the audio level crosses above threshold.
%syllEndTimes: the end time of each syllable.  This is not the offset of
%sound, but the time at which the audio level crosses below threshold.
%noiseEst, noiseStd: the estimated noise level and variance.
%soundEst, soundStd: the estimated signal level and variance.

cd(directory);
d = dir ('*.wav');    
fs=44100;
mkdir('song files');
mkdir('NOT song files');
bDebug =false;
if ~isempty(d)
    for i=1:length(d)
        filename = d(i).name;
        audio = wavread(filename);
        syllStartTimes = [];
        syllEndTimes = [];
        fMinSyllDuration = .01; %secs
        fMinIntervalDuration = .007; %secs
        fMaxSyllDuration = 1.2; %sec

        %utilizes several thresholds and criteria
        %1. syllable trig threshold power - the syllable has to reach this
        %   threshold at least once.
        %2. syllable cont threshold power - all time between this power is
        %   considered the syllable.
        %3. min interval duration.
        %4. min syllable duration.

        %5. max syllable duration

        %get the log power in the audio signal
        [audioLogPow,audio] = getLogPower(audio, fs);

        %Estimate threshold to discrimate between sound and noise using a mixture
        %of two gaussians model.
        [noiseEst, soundEst, noiseStd, soundStd] = aSAP_estimateTwoMeans(audioLogPow);    
        if(noiseEst>soundEst)
            movefile([directory,'\', filename], [directory,'\NOT song files\',filename]);
        else

            %Compute the optimal classifier between the two gaussians...
            p(1) = 1/(2*soundStd^2) - 1/(2*noiseStd^2);
            p(2) = (noiseEst)/(noiseStd^2) - (soundEst)/(soundStd^2);
            p(3) = (soundEst^2)/(2*soundStd^2) - (noiseEst^2)/(2*noiseStd^2) + log(soundStd/noiseStd);
            disc = roots(p);
            disc = disc(find(disc>noiseEst & disc<soundEst));
            if(length(disc)==0)
               movefile([directory,'\', filename], [directory,'\NOT song files\',filename]);
            else
                disc = disc(1);

                %Set the thresholds based on these estimates
                thresSyll = noiseEst + 0.8 * (disc - noiseEst);
                thresTrig = soundEst - 0.6 * (soundEst - disc);
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
                    NFFT=512;
                    syllable=[];
                    for j=1:length(beginSyll)
                        y=fft(audio(beginSyll(j):endSyll(j)),NFFT);
                        power=2*abs(y(1:NFFT/2));
                        f = fs/2*linspace(0,1,NFFT/2);
                        if (bDebug)
                            figure;plot(f,power);
                            soundsc(audio(beginSyll(j):endSyll(j)),fs);
                        end
                        
                        if sum((power(find(f>1000,1,'first'):find(f < 5000,1,'last'))))>sum(power(find(f>5000,1,'first'):find(f < 8000,1,'last')))
                            syllable(j)=1;
                        else
                            syllable(j)=0;
                        end
                        if (bDebug)
                            lowpower=sum(power(find(f>1000,1,'first'):find(f < 5000,1,'last')))
                            highpower=sum(power(find(f>5000,1,'first'):find(f < 8000,1,'last')))
                        end
                    end
                    realSyll=find(syllable==1);
                    beginSyll = beginSyll(realSyll);
                    endSyll = endSyll(realSyll);
                    syllStartTimes = (beginSyll -1)/ fs;
                    syllEndTimes = (endSyll-1) / fs;
                end
%                 if(bDebug)
%                     s2 = subplot(2,1,2);
% %                     [SAP_Feats, m_spec_deriv] = aSAP_generateASAPFeatures(audio, fs);
% %                     aSAP_displaySpectralDerivative(m_spec_deriv, Parameters);
%                     plot([0:length(audio)-1]/fs,audioLogPow);
%                     allCross = [syllStartTimes, syllEndTimes];
%                     for(i = 1:length(allCross))
%                         axes(s1);
%                         line([allCross(i), allCross(i)], ylim, 'Color', 'red');
%                         axes(s2);
%                         line([allCross(i), allCross(i)], ylim, 'Color', 'red');
%                     end
%                     linkaxes([s1,s2], 'x')
%                     figure(gcf);
%                     pause;
%                 end
                num = BoutsFromRawAudio(syllStartTimes, syllEndTimes);
                if (num>0)
                    movefile([directory,'\', filename], [directory,'\song files\',filename]);
                else
                    movefile([directory,'\', filename], [directory,'\NOT song files\',filename]);

                end
            end
        end
    end
end
