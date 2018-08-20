function [newTemp, newsyllBreaks, newmotifBreaks] = platonicWarp()%(template,syllBreaks,motifBreaks,meanPath,pathT)

[fname loc] = uigetfile('C:\Users\Tim\Desktop\For Ali Paper Test Alignments\Pur696\','Select template source');
load([loc fname])
template = template;
syllBreaks = templatesyllBreaks;
motifBreaks = templatemotifBreaks;

[fname loc] = uigetfile('C:\Users\Tim\Desktop\For Ali Paper Test Alignments\Pur696\','Select path source');
load([loc fname])
meanPath = paths.pathMean;
pathT = paths.pathT;

%Create the inverse to the given meanPath
invPath = [pathT (-.5*meanPath'+pathT)];
invPath(invPath(:,2)>size(template,2),2) = size(template,2);
newTemp = alignSeries(template,round(invPath));

figure(228)
% subplot(2,2,1)
% imagesc(-1*template)
% axis xy
% 
% subplot(2,2,3)
% imagesc(-1*newTemp)
% axis xy

subplot(2,2,3:4)
hold on
plot(pathT',meanPath)
plot(pathT',invPath(:,2)'-pathT','r')
hold off

%newsyllBreaks = getRemap(invPath,syllBreaks);
%newmotifBreaks = getRemap(invPath,motifBreaks);
[newsyllBreaks,newmotifBreaks] = getEdges(newTemp,syllBreaks,[],[]);

subplot(2,2,1:2)
imagesc(-1*template)
hold on
for i = 1:length(newsyllBreaks(:))
    line([newsyllBreaks(i), newsyllBreaks(i)],[0,size(template,1)],'Color','k','LineWidth',2)
end
axis tight; axis xy;
hold off




function [warpedOut] = getRemap(path,anchors)
[m, n] = size(anchors);
warpedOut = zeros(m,n);

for i = 1:m
    for j = 1:n
        ind = find(path(:,1)==anchors(i,j));
        %ind = find(path(:,1)>=(anchors(i,j)-43) & path(:,1)<=anchors(i,j));
        warpedOut(i,j) = round(mean(path(ind,2)));
    end
end

function [syllBreaks,motifBreak] = getEdges(data,starts,buffer,type)

Crossings.Up = [];
Crossings.Down = [];
    
%Sum across frequency bins for each rendition
specPowTS = -1*sum(data);
specPow = sum(sum(data));

%Rerun EM until a solution is found
iterations = 1;
gMix = [];

while (isempty(Crossings.Up) || isempty(Crossings.Down))
    gMix = gmdistribution.fit(specPowTS',2);
    [powerM,pnt] = sort(gMix.mu);
    temp = sqrt(squeeze(gMix.Sigma));
    powerStd = temp(pnt);
    thresh = powerM(1)+2*powerStd(1); %1SDs above the silent Gaussian mean

    %Find the crossing points to estimate the onset and offsets of
    %motifs
    ind = 2:length(specPowTS);
    Crossings.Up = find(specPowTS(ind)>thresh & specPowTS(ind-1)<thresh);
    Crossings.Down = find(specPowTS(ind)<thresh & specPowTS(ind-1)>thresh);
    
    if iterations >5
        print(['gmdistribution.fit has been run ' double2str(iterations) ' without converging on a solution.'])
    end
    iterations = iterations + 1;
end

%Use the starting positions to capture the start/stop points for all
%syllables of the rendition
syllBreaks = [];
motifBreak = [];

onStarts = starts(:,1);
offStarts = starts(:,2);

 for k = 1:length(onStarts)
     [offset, pntr] = min(abs(Crossings.Up-onStarts(k)));
      if offset > 10 && strcmp(type,'robust')
         syllBreaks(k,1) = onStarts(k);
      else
        syllBreaks(k,1) = Crossings.Up(pntr);
      end
     
     [offset, pntr] = min(abs(Crossings.Down-offStarts(k)));
      if offset > 10 && strcmp(type,'robust')
         syllBreaks(k,2) = offStarts(k);
      else
        syllBreaks(k,2) = Crossings.Down(pntr);
      end
 end
motifBreak(1) = syllBreaks(1,1); motifBreak(2) = syllBreaks(end,2);