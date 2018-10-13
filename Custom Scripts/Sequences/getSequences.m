function [allSeq, startSeq, counter]=getSequences(elements, numberSyll, startends)
%Gets all the sequences containing numSyll.  Note that this is entirely
%dependent on the current filter settings.

%Initialize placeholders
song_interruptions=[];
allSeq=[];
startSeq=[];
counter=1;

%Restrict sequence search to those files filtered by the user (legacy)
filt_elements = elements;

%Cycle through records to build up a map of the sequences in the annotation
for i = 1:length(filt_elements)
    
    %Find all pauses greater than 200ms to determine boundaries of song motifs
    for j = 1:length(filt_elements{i}.segFileStartTimes)-1
        %Calculate the separation between adjacent sound segments
        pauses(j) = filt_elements{i}.segFileStartTimes(j+1)-filt_elements{i}.segFileEndTimes(j);
        
        %If that seaparation is greater than the threshold (200ms), count
        %it as a genuine break in singing
        song_interruptions=find(pauses>0.2);
    end
    
    %If you're running only 2-syl sequences and startends are requested,
    %insert the flags (-98=start; -99=end)
    if startends && numberSyll==2
        ins = [-99, -98];
        segTypes = -98;
        for j = 1:numel(filt_elements{i}.segType)
            %In to the list
            segTypes = [segTypes, filt_elements{i}.segType(j)];
            
            %append a flag
            if ismember(j, song_interruptions)
                segTypes = [segTypes, ins];
            end
        end
        segTypes = [segTypes, -99];
    else
        segTypes = filt_elements{i}.segType;
    end
    
    
    %If there are enough syllables in the file, start it up
    if length(segTypes)>=numberSyll
        
        %Step through the annotated intervals and extract the subsequent sylTypes
        %for j = 1:length(filt_elements{i}.segFileStartTimes)-numberSyll+1
        for j = 1:length(segTypes)-numberSyll+1
            %Grab the current sequence
            %current_seq(1:numberSyll) = filt_elements{i}.segType(j:j+numberSyll-1);
            current_seq(1:numberSyll) = segTypes(j:j+numberSyll-1);
            
            
            
            %If the sequence has valid labels and doesn't include the
            %previously identified song breaks, put it in the stack
            if isempty(find(current_seq==-1 | current_seq>9)) && startends && numberSyll==2 && ~isequal(current_seq, ins)
                %Throw it on the stack
                allSeq(counter,1:numberSyll)=current_seq;
                
                %Note the sequence origin (annotation file index, position of lead syllable in annotation)
                startSeq(counter,:)=[i,j];
                
                %Update counter
                counter=counter+1;
                
                %There is a song break in this segment, so add that into the running stack
            elseif (isempty(find(current_seq<1 | current_seq>9)) && ~startends && isempty(find(song_interruptions>j-1 & song_interruptions<j+numberSyll-1)))
                %Throw it on the stack
                allSeq(counter,1:numberSyll)=current_seq;
                
                %Note the sequence origin (annotation file index, position of lead syllable in annotation)
                startSeq(counter,:)=[i,j];
                
                %Update counter
                counter=counter+1;
            end
        end
    end
end
