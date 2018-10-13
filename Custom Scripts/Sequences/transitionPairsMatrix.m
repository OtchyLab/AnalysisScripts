function seqs = transitionPairsMatrix(annots, startends, plotIt)

%Defaults
if nargin < 2
    startends = false;
    plotIt = false;
end

if nargin < 3
    plotIt = false;
end
    
%Loop through the list of files and copy data from the harddrive to the
%local variable
numSyl=2; %sequence depth
seqs = [];
for i = 1:numel(annots)
    %Select the elements node
    t = annots{i};
    
    %Retrieve sequences
    [allSeq, ~, counter]=getSequences(t, numSyl, startends);
    
    %Determine what syllables types are present
    sylTypes = sort(unique(allSeq(:)));
    realSyls = sylTypes(sylTypes>0 & sylTypes <10);
    axesLabels = realSyls;
    
    %Build mapping
    if ismember(-98, sylTypes)
        map = [-98; realSyls; -99];
        axesLabels = {'S', axesLabels, 'E'};
    end
    
    %Shuffle into a count matrix (rows=from; columns=to)
    seqMatrix=zeros(numel(sylTypes), numel(sylTypes));
    for j=1:counter-1
        from = find(map==allSeq(j,1));
        to = find(map==allSeq(j,2));
        seqMatrix(from,to)=seqMatrix(from,to)+1;
    end
    
    %Copy out to the carried variable
    seqs(i,:,:) = seqMatrix;
    
    if plotIt
        %Plot the figure
        figure;
        h=imagesc([1 numel(sylTypes)],[1 numel(sylTypes)], seqMatrix);
        xlabel('To'); ylabel('From')
        set(gca, 'Box', 'off', 'TickDir', 'out')
        set(gca, 'XTick', 1:numel(sylTypes), 'XTickLabels', axesLabels, 'YTick', 1:numel(sylTypes), 'YTickLabels', axesLabels)
    end
end










