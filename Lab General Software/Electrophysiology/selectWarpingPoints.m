function cell = selectWarpingPoints(cell, templateAudio, templateStartTime, templateWarpMarks)

limitBuffer = [-.1,+.1];
sampRate = cell.exper.desiredInSampRate;

%Plot the template...
h = figure(1006);
subplot(4,1,1)
displayAudioSpecgram(templateAudio, sampRate, templateStartTime);
axis tight;
xlim(xlim + limitBuffer);
limits = xlim;

subplot(4,1,2);
logpow = log(templateAudio.^2);
time = ((0:length(templateAudio)-1)/sampRate) + templateStartTime;
plot(time, logpow);
axis tight;
xlim(limits)

for(nWarp = 1:length(templateWarpMarks))
    markTime = ndx2time(templateWarpMarks(nWarp)-1, sampRate) + templateStartTime;
    subplot(4,1,1); hold on; l = line([markTime,markTime], ylim); set(l,'Color','red');
    subplot(4,1,2); hold on; l = line([markTime,markTime], ylim); set(l,'Color','red');   
end

%Compute normalization factor
norm = prctile(templateAudio,95);

for(nMotif = 1:length(cell.motifData))
    %Normalize this motif to the template.
    
    motifData = cell.motifData{nMotif};
    if(~isfield(motifData,'warpMarks'))
        audio = motifData.audio;
        scale = norm/prctile(audio,95);
        normaudio = audio*scale;
		startTime = ndx2time(motifData.marker,sampRate);    
        
        figure(h);
        subplot(4,1,3)
        displayAudioSpecgram(normaudio, sampRate, -startTime);
        axis tight;
        xlim(limits);
        title(['Motif #',num2str(nMotif)]);
        
        subplot(4,1,4);
        logpow = log(normaudio.^2);
        time = ((0:length(normaudio)-1)/sampRate) - startTime;
        plot(time, logpow);
        axis tight;
        xlim(limits);
    	title(['Motif #',num2str(nMotif)]);
    
        threshold  = -11;
        warpMarks = [];
        warpHandles = [];
        count = 0;
        while(true)
           pause(.05);
           [x,y,opt] = ginput(1)
           if(opt == '.')
               cell.motifData{nMotif}.warpMarks = sort(warpMarks);
               cell.motifData{nMotif}.oldMarker = motifData.marker;
               cell.motifData{nMotif}.marker = warpMarks(3);
               break;
           elseif(opt == 3)
                %Find marker location in motif.
                count = count + 1;
                aproxMark = max(round((x+startTime) * sampRate) + 1,1);
                mark = find(logpow(aproxMark:end) > threshold);
                mark = mark(1) + aproxMark-1;
                warpMarks(count) = mark;
                markTime = ndx2time(mark-motifData.marker-1, sampRate);
                subplot(4,1,3); warpHandles(count,1) = line([markTime,markTime], ylim); set(warpHandles(count,1),'Color','red');
                subplot(4,1,4); warpHandles(count,2) = line([markTime,markTime], ylim); set(warpHandles(count,2),'Color','red');
            elseif(opt == 1)
                count = count + 1;
                aproxMark = min(round((x+startTime) * sampRate) + 1,length(logpow));
                mark = find(logpow(1:aproxMark) > threshold);
                mark = mark(end);
                warpMarks(count) = mark;
                markTime = ndx2time(mark-motifData.marker-1, sampRate);
                subplot(4,1,3); warpHandles(count,1) = line([markTime,markTime], ylim); set(warpHandles(count,1),'Color','red');
                subplot(4,1,4); warpHandles(count,2) = line([markTime,markTime], ylim); set(warpHandles(count,2),'Color','red');     
            elseif(opt == 'd')
                warpMarks = warpMarks(1:end-1);
                delete(warpHandles(count,:));
                warpHandles = warpHandles(1:end-1,:);
                count = count -1;
            elseif(opt == 't')
                threshold = input('Enter a threshold?');
            elseif(opt == 's');
                saveCell(cell);
            end
        end        
    end    
    if(mod(nMotif,20) == 0)
        saveCell(cell);
    end
end

saveCell(cell);

