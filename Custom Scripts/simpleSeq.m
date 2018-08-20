function [sylStream, timestamp,allSeq,allSeqTime] = simpleSeq(elements)
%go through the elements structure of an annotation file and whenever there is a stretch of contiguous
%sequence assign the transition to the matrix.

%Reset placeholders
counter=1;
numberSyll=2;
boutbreak = 99;
sylStream = [];
timestamp = [];
%Find all pauses greater than 100ms to determine boundaries of song motifs
for i = 1:length(elements)
    pauses = [];
    for j = 1:length(elements{i}.segFileStartTimes)-1
        pauses(j) = elements{i}.segFileStartTimes(j+1)-elements{i}.segFileEndTimes(j);
        song_interruptions=find(pauses>0.1);
    end

    if length(elements{i}.segFileStartTimes)>=numberSyll
        %Reshape the syllable sequence for the whole file to reflect the
        %start and ends of song bouts
        segTypeExt = elements{i}.segType;
        segTime = elements{i}.segAbsStartTimes';
        for j = 1:length(song_interruptions)
%             ind = song_interruptions(j)+(2*(j-1));
            ind = song_interruptions(j)+((j-1));
            
            pre = segTypeExt(1:ind);
            post = segTypeExt((ind+1):end);
            segTypeExt = [pre; boutbreak; post];
            
            preT = segTime(1:ind);
            postT = segTime((ind+1):end);
            segTime = [preT; preT(end); postT];
        end
        segTypeExt = [segTypeExt; 99];
        segTime = [segTime; segTime(end)];
        
        for j = 1:length(segTypeExt)-numberSyll+1
            current_seq(1:numberSyll) = segTypeExt(j:j+numberSyll-1);
%             if (isempty(find(current_seq<1 | current_seq>9)) && isempty(find(song_interruptions>j-1 & song_interruptions<j+numberSyll-1)))
%             if (isempty(find(current_seq==-1 | current_seq>99)) && ~isequal(current_seq,boutbreak'))
%             if (~isequal(current_seq,boutbreak'))
%             if isempty(find(song_interruptions>j-1 & song_interruptions<j+numberSyll-1)) %Only exclude song breaks... all other transitions included
                allSeq(counter,1:numberSyll)=current_seq;
                allSeqTime(counter) = segTime(j);
                startSeq(counter,:)=[i,j];
                counter=counter+1;
            
        end
        
        sylStream = [sylStream; segTypeExt];
        timestamp = [timestamp;segTime];
    end
end
