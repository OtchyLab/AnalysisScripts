%function sequenceStats_annot
%Function will calculate and plot a few different representations of
%sequence variability. 
%
%The main input folder is selected by the user using a GUI prompt. in that
%folder should be several saved annotation files created by the Talon
%applications (i.e., TweetVision or TweetVisionLite). Each file in that 
%folder should be named:
%     BirdName_YYMMDD_annotation.mat
%
% Created 10/02/2018 TMO

%% Set Defaults
startends = true; %Include start and end points in sequences
plotIt = true; %Boolean flag to start plot the output

surgeryDate = '180111'; %change this as necessary

%% Load from file

%Select the data folder to process
parent = uigetdir(cd, 'Which folder to process?');

%Get the list of files to process
annotFiles = dir([parent, filesep, '*annotation.mat']);

%Figure out how many files you're going to process
numFiles = numel(annotFiles);

%Loop through the list of files and copy data from the harddrive to the
%local variable
annots = [];
for i = 1:numFiles
    %Select file to load
    filename = annotFiles(i).name;

    %Load from the file
    load([parent, filesep, filename], 'elements');
    
    %Copy to local variable
    annots{i} = elements;
    
    %Delete the loaded variable
    clear('elements')
end

%% Extract date information for later plotting
surgeryNum = datenum(surgeryDate, 'yymmdd');
days = [];
for i = 1:numFiles
    %Select file to load
    filename = annotFiles(i).name;

    %Split up the filename using the underscore
    sp = regexp(filename, '_', 'split');
    
    %Select the portion that includes the date
    dStr = sp{2};
    
    %Convert to a datenum for ease
    dNum = datenum(dStr, 'yymmdd');
    
    %Day before or since surgery
    days(i) = dNum - surgeryNum;

end

%% Run the sequence extraction steps

%Create sequence pairs matrices for each daya
%seqPairs = transitionPairsMatrix(annots, startends, plotIt);

%Get the distibution of sequences
seqTokens= [];
seqCounts = [];
rankRange = 1:8;
for i = 1:numFiles
    for j = rankRange
        %Capture all squences of rank j
        [allSeq, startSeq] = sortSequences(annots{i}, j);
        
        %Sort sequences by prevalance
        histx=(1:10^j);
        hist_seq=hist(allSeq,histx);
        [sortSeq, sortIndex]=sort(hist_seq,'descend');
        numSeq=length(find(sortSeq>0));
        
        %Pad if short
        if numSeq<100
            tsS = zeros(1,100);
            tsS((1:numSeq)) = sortSeq(1:numSeq);
            sortSeq = tsS;
            
            tsI = zeros(1,100);
            tsI((1:numSeq)) = sortIndex(1:numSeq);
            sortIndex = tsI;
        end
        
        %Dump to a carried variable
        seqTokens(i,j,:) = sortIndex(1:100);
        seqCounts(i,j,:) = sortSeq(1:100);
        
    end
    a = 1;
end

