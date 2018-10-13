function [allSeqNumbers, startSeq]=sortSequences(elements, numberSyll)
%make a cell array with all the song sequences

[allSeq, startSeq, counter]=getSequences(elements, numberSyll, false);
if (allSeq==0)
    warndlg('There are no sequences in the specified range.');
    uiwait;
    return;
end
%go through this cell array and build up an array of sequences
%numberSyll long
for i=1:numberSyll
    multiplier(i)=10^(numberSyll-i);
end
allSeqNumbers=[];
for i=1:counter-1
    allSeqNumbers(i)=dot(allSeq(i,:),multiplier);
end