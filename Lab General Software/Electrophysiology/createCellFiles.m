function cell = createCellFiles(cell, fileNumbers)

load([cell.dir,cell.filename]);
sampRate = cell.exper.desiredInSampRate;

if(~isfield(cell,'fileNumbers'))
    cell.fileNumbers = [];
end
cell.fileNumbers = union(cell.fileNumbers, fileNumbers);

%If new cell initialize.
if(~isfield(cell,'stimFileNums'))
    cell.stimFileNums = [];
end

if(~isfield(cell,'motifData'))
    numMotifData = 0;
else
    numMotifData = length(cell.motifData);
end

if(~isfield(cell,'messyData'))
    numMessyData = 0;
else
    numMessyData = length(cell.messyData);
end

if(~isfield(cell,'spontData'))
    numSpontData = 0;
else
    numSpontData = length(cell.spontData);
end

for(fileNum = fileNumbers) 
    %Load the data to be processed.
    audio = loadAudio(cell.exper, fileNum);
    audio = audio - mean(audio); %demean audio
    sig = loadData(cell.exper, fileNum, cell.chanNum);
    time = [1:length(audio)] / sampRate;
    
    %Display the data
    figure(1001);
    subplot(3,1,1);
    displayAudioSpecgram(audio, sampRate); axis tight;
    title(['file',num2str(fileNum),' chan',num2str(cell.chanNum)]);
    subplot(3,1,2);
    plot(time, audio);  axis tight;
    subplot(3,1,3);
    plot(time, sig);  axis tight;
    
    while(true)        
        %Wait for feedback
        leftButton = 1;
        rightButton = 2;
        
        [timeStart, y, opt] = ginput(1);
        if(opt == leftButton)
            [timeStop,y, opt] = ginput(1);
            ndxStart = floor(timeStart*sampRate);
            ndxStop = floor(timeStop*sampRate);
            
            %compute log power of audio
            logpow = log(audio(ndxStart:ndxStop).^2);
            
            %Display clipped out version.
            figure(1002);
            subplot(4,1,1);
            displayAudioSpecgram(audio(ndxStart:ndxStop), sampRate, ndxStart/sampRate); axis tight;
            subplot(4,1,2);
            plot(time(ndxStart:ndxStop), audio(ndxStart:ndxStop)); axis tight;
            subplot(4,1,3);
            plot(time(ndxStart:ndxStop), logpow); axis tight;
            subplot(4,1,4);
            plot(time(ndxStart:ndxStop), sig(ndxStart:ndxStop));  axis tight;
            
            %Find marker location in motif.
            [aproxMarkerTime, y, opt] = ginput(1);
            aproxMarker = round(aproxMarkerTime * sampRate);
            marker = find(logpow(aproxMarker - ndxStart:end) > -10);
            marker = marker(1) + aproxMarker - ndxStart -1;
            for(i = 2:4)
                figure(1002);
                h = subplot(4,1,i);
                yLimits = get(h,'yLim');
                markTime = (marker + ndxStart) / sampRate;
                h = line([markTime, markTime],yLimits)
                set(h,'Color','red');
            end
            
            %Let user zoom, play around ect... until y or n is pressed.            
            h = figure(1002);
            title('Is this good data? (y or n)?');
            set(h,'CurrentCharacter',' ');
            while(true)
                pause(.05);
                charGoodData = get(h,'CurrentCharacter');
                if(charGoodData == 'y' | charGoodData == 'n')
                    break;
                end
            end
            if(charGoodData == 'y')
                numMotifData = numMotifData + 1;
                cell.motifData{numMotifData}.audio = audio(ndxStart:ndxStop);
                cell.motifData{numMotifData}.sig = sig(ndxStart:ndxStop);
                cell.motifData{numMotifData}.origFile = fileNum;
                cell.motifData{numMotifData}.ndxStart = ndxStart; %Relative to file start.
                cell.motifData{numMotifData}.ndxStop = ndxStop; %Relative to file start.
                cell.motifData{numMotifData}.marker = marker; %Relative to index start.
                
                %Display black line on main figure to mark parsed region.
                figure(1001); subplot(3,1,1);
                hold on;
                line([ndxStart/sampRate, ndxStop/sampRate], [2000,2000]);
                hold off;
            else
                numMessyData = numMessyData + 1;
                cell.messyData{numMessyData}.audio = audio(ndxStart:ndxStop);
                cell.messyData{numMessyData}.sig = sig(ndxStart:ndxStop);
                cell.messyData{numMessyData}.origFile = fileNum;
                cell.messyData{numMessyData}.ndxStart = ndxStart;
                cell.messyData{numMessyData}.ndxStop = ndxStop;
                cell.messyData{numMessyData}.marker = marker; 
                
                %Display black line on main figure to mark parsed region.
                figure(1001); subplot(3,1,1);
                hold on;
                line([ndxStart/sampRate, ndxStop/sampRate], [2000,2000]);
                hold off;
            end
            
            figure(1001);
        elseif(opt == rightButton) %Spontaneous activity.
            [timeStop,y, opt] = ginput(1);
            ndxStart = floor(timeStart*sampRate);
            ndxStop = floor(timeStop*sampRate);
            numSpontData = numSpontData + 1;
            cell.spontData{numSpontData}.audio = audio(ndxStart:ndxStop);
            cell.spontData{numSpontData}.sig = sig(ndxStart:ndxStop);
            cell.spontData{numSpontData}.origFile = fileNum;
            cell.spontData{numSpontData}.ndxStart = ndxStart;
            cell.spontData{numSpontData}.ndxStop = ndxStop;
            cell.spontData{numSpontData}.marker = marker - ndxStart;    
        elseif(opt == 's') %file has stims
            %Let user zoom, play around ect... until spacebar is pressed
            h = figure(1001);
            zoom on;
            set(h,'CurrentCharacter','s');
            while(true)
                pause(.05);
                char = get(h,'CurrentCharacter');
                if(char == ' ')
                    break;
                end
            end
            figure(1001);
            %Then save stim filenum in cell structure.
            cell.stimFileNums = [cell.stimFileNums, fileNum];
        elseif(opt == '.') %next file
            break;   
        elseif(opt == 'c')
            if(~isfield(cell,'comments'))
                nComment = 1;
            else
                nComment = length(cell.comments) + 1;
            end
            cell.comments{nComment} = input('Comment: ','s');
        elseif(opt == 'q') %next file
            return;  
        elseif(opt == 'z')
            %Let user zoom, play around ect... until spacebar is pressed.
            h = figure(1001);
            zoom on;
            set(h,'CurrentCharacter','z');
            while(true)
                pause(.05);
                char = get(h,'CurrentCharacter');
                if(char == ' ')
                    break;
                end
            end
            figure(1001);
        end
    end
    saveCell(cell);
end

cell.creationComments = input('Comments on the cell?', 's');
saveCell(cell);