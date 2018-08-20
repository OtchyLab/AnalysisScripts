function cell = extractNonSongSpikeTrains(cell)

cell = loadCell(cell.num);

for nNonSong = 1:length(cell.nonSongData)
    spikeNdx = [];
    spikeMarks = []
    nonSongData = cell.nonSongData{nNonSong};
    
    leftButton = 1;
    rightButton = 3;
    if(~isfield(nonSongData, 'spikeNdx'))
        h = figure(1003);
        clf;
        plot(nonSongData.sig);
        title(['Cell:', num2str(cell.num), ' Motif:', num2str(nNonSong), 'of ', num2str(length(cell.nonSongData))]);
        hold on;
        [x,threshold, opt] = ginput(1);
        pause(.01);
        
        if(opt == leftButton | opt == rightButton)
            sig = nonSongData.sig;
            
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
                    warning(['Spikewidth: ', num2str(spikeWidth), ' seems to high.']);
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
                    spikeMarks(nSpike) = plot(spikeNdx(nSpike), peakV, 'ro');
	
                end
            else
                spikeNdx = [];
            end
            
            figure(1002);
            %Wait for quick sanity check approval.
            h = figure(1003);
            set(h,'CurrentCharacter',' ');
            while(true);                pause(.05);
                char = get(h,'CurrentCharacter');
                set(h,'CurrentCharacter','~');
                if(char == 'y')
                    cell.nonSongData{nNonSong}.spikeNdx = spikeNdx;
                    cell.nonSongData{nNonSong}.bSpikesHandChecked = false;
                    cell.nonSongData{nNonSong}.spikeMethod = 'threshold';  
                    cell.nonSongData{nNonSong}.spikeThreshold = threshold;
                    break;     
                elseif(char == 'a')
                    figure(h);
                    [x,y,opt] = ginput(1);
                    x = round(x)
                    [peak, spike] = max(nonSongData.sig(x-200:x+200));
                    spike = spike + x - 50 -1;
                    spikeNdx = sort([spikeNdx,spike]);
                    location = find(spikeNdx == spike);
                    mark = plot(spike, peak, 'ro');
                    spikeMarks = [spikeMarks(1:location-1),mark,spikeMarks(location:end)];
                    figure(h);
                    pause(.1);
                elseif(char == 'd');
                    figure(h);
                    [xStart,yStart,opt] = ginput(1);
                    [xStop,yStop,opt] = ginput(1);
                    delSpikeNums = find(spikeNdx > xStart & spikeNdx < xStop);
                    keepSpikeNums = find(spikeNdx <= xStart | spikeNdx >= xStop);
                    delete(spikeMarks(delSpikeNums));
                    spikeMarks = spikeMarks(keepSpikeNums); 
                    spikeNdx = spikeNdx(keepSpikeNums);  
                    figure(h);
                    pause(.1);
                elseif(char == 'n')
                    break;
                end
            end                 
        elseif(opt == '.' | opt== 'n')
            %do nothing and go to next file.
        end
        clf(1003);
    end
end

saveCell(cell);
        
                
            
         
        