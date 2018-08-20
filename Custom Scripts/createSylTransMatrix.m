function transMatrix = createSylTransMatrix(Elements)
%Elements is simply the 'elements' data from the annotation files. No need in load any actual song recordings.

%Generate sequence transition probabilities
[transMatrix, sylTypes] = SequenceMatrixSE(Elements);


%Display those
h = figure;
axes1 = axes('Parent',h);%'YDir','reverse','XTick',[1 2 3 4 5],'Layer','top');
imagesc(transMatrix);
hold on
axis tight

%Plot matrix
for i = 1:length(sylTypes)
    tempTypes{i} = num2str(sylTypes(i));
    totFol = sum(transMatrix(i,:));
     for j = 1:length(sylTypes)
         text(j-.1,i,num2str(transMatrix(i,j)/totFol,2));
     end
end
tempTypes{sylTypes == -99} = 'S';
tempTypes{sylTypes == 99} = 'E';
set(axes1,'XTickLabel',tempTypes,'XTick',1:length(sylTypes))
set(axes1,'YTickLabel',tempTypes,'YTick',1:length(sylTypes))
hold off


function [transMatrix,sylTypes] = SequenceMatrixSE(elements)
%go through the annotations and whenever there is a stretch of contiguous
%sequence assign the transition to the matrix.

%Reset placeholders
counter=1;
song_interruptions=[];
numberSyll=2;
boutbreak = [99;-99];

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
        for j = 1:length(song_interruptions)
            ind = song_interruptions(j)+(2*(j-1));
            pre = segTypeExt(1:ind);
            post = segTypeExt((ind+1):end);
            segTypeExt = [pre; boutbreak; post];
        end
        segTypeExt = [-99; segTypeExt; 99];
        
        for j = 1:length(segTypeExt)-numberSyll+1
            current_seq(1:numberSyll) = segTypeExt(j:j+numberSyll-1);
            %if (isempty(find(current_seq<1 | current_seq>9)) && isempty(find(song_interruptions>j-1 & song_interruptions<j+numberSyll-1)))
%             if (isempty(find(current_seq==-1 | current_seq>99)) && ~isequal(current_seq,boutbreak'))
            if (~isequal(current_seq,boutbreak'))
                allSeq(counter,1:numberSyll)=current_seq;
                startSeq(counter,:)=[i,j];
                counter=counter+1;
            end
        end
    end
end

%numberSyllTypes=max(max(allSeq));
numberSyllTypes=length(unique(allSeq));
sylTypes = sort(unique(allSeq));
transMatrix=zeros(numberSyllTypes,numberSyllTypes);   
for i=1:counter-1
    from=find(sylTypes==allSeq(i,1));
    to=find(sylTypes==allSeq(i,2));
    transMatrix(from,to)=transMatrix(from,to)+1;
end