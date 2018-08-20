function cell = createNonSongCellFiles(cell)

load([cell.dir,cell.filename]);
sampRate = cell.exper.desiredInSampRate;

if(~isfield(cell,'nonSongData'))
    numNonSongData = 0;
else
    numNonSongData = length(cell.nonSongData);
end

%Because first few files are typically bad... skip them.
fileNumbers = cell.fileNumbers;
%if(length(fileNumbers > 15))
%    fileNumbers = fileNumbers(10:end);
%end

fileNumbers = setdiff(fileNumbers, cell.stimFileNums);
subSamp = 10;
for(fileNum = fileNumbers) 
    %Load the data to be processed.
    audio = loadAudio(cell.exper, fileNum);
    redAudio = audio(1:subSamp:end);
    audio = audio - mean(audio); %demean audio
    time = [1:length(audio)] / sampRate;
    sig = [];
    
    %Display the data
    figure(1001);
    subplot(1,1,1);
    plot(time(1:subSamp:end), redAudio);  axis tight;
    title(['cell',num2str(cell.num),' file',num2str(fileNum),' chan',num2str(cell.chanNum)]);
    
    %Wait for feedback
    leftButton = 1;
    rightButton = 2;    
    
    while(true)           
        [timeStart, y, opt] = ginput(1);
        if(opt == leftButton)
            [timeStop,y, opt] = ginput(1);
            ndxStart = floor(timeStart*sampRate);
            ndxStop = floor(timeStop*sampRate);
            beep;

            if(length(sig) == 0)
                sig = loadData(cell.exper, fileNum, cell.chanNum);
            end
            
            numNonSongData = numNonSongData + 1;
            cell.nonSongData{numNonSongData}.audio = audio(ndxStart:ndxStop);
            cell.nonSongData{numNonSongData}.sig = sig(ndxStart:ndxStop);
            cell.nonSongData{numNonSongData}.origFile = fileNum;
            cell.nonSongData{numNonSongData}.ndxStart = ndxStart;
            cell.nonSongData{numNonSongData}.ndxStop = ndxStop;             
            
            %Display black line on main figure to mark parsed region.
            figure(1001); subplot(1,1,1);
            hold on;
            line([ndxStart/sampRate, ndxStop/sampRate], [.02,.02]);
            hold off;
        elseif(opt == '.') %next file
            break;   
        elseif(opt == 'q')
            return;  
        end
    end
    saveCell(cell);
end