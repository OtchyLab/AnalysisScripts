function spikeNdx=extractSpikeTrainsBence(signal)

    
    spikeNdx = [];
    motifData.sig = signal;
    
    leftButton = 1;
    rightButton = 3;
%    if(~isfield(motifData, 'spikeNdx'))
        figure(1003);
        plot(motifData.sig);
        hold on;
        [x,threshold, opt] = ginput(1);
        pause(.01);
        
        if(opt == leftButton | opt == rightButton)
            sig = motifData.sig;
            
            %Get Threshold Crossings
            if(opt == leftButton) %positive threshold;
                aboveThres = find(sig > threshold);
            else %negative threshold.
                aboveThres = find(sig < threshold);
            end
            
            if(length(aboveThres == 0))
                downCrossNdx = find(diff(aboveThres) > 1);
                upCrossNdx = downCrossNdx + 1;
                downCross = [aboveThres(downCrossNdx); aboveThres(end)];
                upCross = [aboveThres(1); aboveThres(upCrossNdx)];
                
                %Find spike width for faster max finding...
                spikeWidth = max(downCross - upCross);
                if(spikeWidth > 15)
                   % warning(['Spikewidth: ', num2str(spikeWidth), ' seems to high.']);
                end
                
                %Find spike peaks, could do without loop if you assume all
                %spikes are approximately the same width, but decided to 
                %code it the slow way... more robust.
                for(nSpike = 1:length(upCross))
                    if(opt == leftButton) %positive threshold;
                        [peakV, peakNdx] = max(sig(upCross(nSpike):downCross(nSpike)));
                        spikeNdx(nSpike) = upCross(nSpike) + peakNdx - 1;
                    else %negative threshold.
                        [peakV, peakNdx] = min(sig(upCross(nSpike):downCross(nSpike)));
                        spikeNdx(nSpike) = upCross(nSpike) + peakNdx - 1;
                    end    
                    figure(1003);
                    plot(spikeNdx(nSpike), peakV, 'ro');
	
                end
            else
                spikeNdx = [];
            end
            
            figure(1002);
            %Wait for quick sanity check approval.
            h = figure(1003);
            set(h,'CurrentCharacter',' ');
            
                           
        elseif(opt == '.' | opt== 'n')
            %do nothing and go to next file.
        end
        %clf(1003);
   %end
end


        
                
            
         
        