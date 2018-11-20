function [seqMatrix, seqProbMatrix, types, mainTrans, h] = makeTransMatrix(annotation)
%Simple function to load an annotation (from passed file location), and
%calculate a transition probability matrix. Outputs are a (1) passed out
%and (2) plotted transition matrix.

%Check the file location
if ~exist(annotation, 'file')
    disp('File wasn''t found... try again')
    return
end

%Load from the file
load(annotation, 'elements');

%Retrieve the list of 2-syllable sequences
[allSeq, ~, counter]=getSequences(elements, 2, false);

%ID unique syllable tokens
types = sort(unique(allSeq));

%Create the transition counts matrix
numberSyllTypes=max(max(allSeq));
seqMatrix=zeros(numberSyllTypes,numberSyllTypes);   
for i=1:counter-1
    from=allSeq(i,1);
    to=allSeq(i,2);
    seqMatrix(from,to)=seqMatrix(from,to)+1;
end

%Count instances of "Froms"
cnts = sum(seqMatrix,2);

%Oops! there are empty rows
if min(cnts) == 0
    %row ranks
    mask = cnts~=0;
    
    %Empty rows remove
    seqMatrix = seqMatrix(mask, mask);
    numberSyllTypes = size(seqMatrix,1);
end

%Recount instances of "Froms"
cnts = sum(seqMatrix,2);

%Calculate sequence probability
seqProbMatrix = seqMatrix./cnts;

%Figure out the main transitions
[m, I] = max(seqProbMatrix,[],2);
mainTrans = [(1:numberSyllTypes)', I];

res = abs(seqProbMatrix-m);
I = find(res>0 & res<0.1);
[row, col] = ind2sub(size(seqProbMatrix), I);
mainTrans = [mainTrans; row, col];

%Prep the figure
h = figure('Name', annotation(1:end-15)); clf
set(gcf, 'Units', 'Inches', 'Position', [1.5, 6, 11, 4.5])

%Plot the counts matrix
subplot(1,2,1); cla
imagesc([1 numberSyllTypes],[1 numberSyllTypes], seqMatrix);
colormap(bone)
colorbar

axis square
xlabel('To')
ylabel('From')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', 1:numel(types), 'XTickLabel', types)
set(gca, 'YTick', 1:numel(types), 'YTickLabel', types)

%Plot the counts matrix
subplot(1,2,2); cla
imagesc([1 numberSyllTypes],[1 numberSyllTypes], seqProbMatrix, [0, 1]);
colormap(bone)
colorbar

axis square
xlabel('To')
ylabel('From')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', 1:numel(types), 'XTickLabel', types)
set(gca, 'YTick', 1:numel(types), 'YTickLabel', types)






