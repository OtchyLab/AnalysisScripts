function displayCellIndexPage(cell)

%Get Motif Counts
numMotifs = 0;
numGoodMotifs = 0;
if(isfield(cell,'motifData'))
    numMotifs = length(cell.motifData);
    for(nMotif = 1:numMotifs)
        if(isfield(cell.motifData{nMotif},'spikeNdx'))
            numGoodMotifs = numGoodMotifs+1;
        end
    end
end

%Prepare the text:
lines1 = sprintf('cell#%g         labbook#%g        bird:%s          exper:%s          signal#:%s',cell.num, cell.labNotebookNum, cell.exper.birdname, cell.exper.expername, cell.signalNum);
lines2 = sprintf('AP:%g LM:%g DV:%g          NumMotifs:%g/%g',cell.AP,cell.LM,cell.DV,numGoodMotifs,numMotifs);
lines3 = sprintf('Identified??:  %g     Current:%g         Est. Latency:%g',cell.antidromicIdentified, cell.antidromicAproxCurr, cell.antidromicGuessLat);
lines4 = cell.antidromicComments;
if(isfield(cell,'creationComments'))
    lines5 = cell.creationComments;
else
    lines5 = '';
end
lines6 = '';
if(isfield(cell,'comments'))
    for(i=1:length(cell.comments))
        lines6 = [lines6,'/',cell.comments{i}];
    end
end
 
numVert = 10;
numHorz = 2;

%Create our figure
h = figure;
set(h,'Renderer','Painters');

%Plot the text
s1 = subplot(numHorz, numVert,[1:4])
u = uicontrol;
set(u,'style','text');
set(u,'units','normalized');
set(u,'horizontalalignment','left');
set(u,'position',[.02,1-2/numVert+.02,.96,2/numVert-.04]);
set(u,'fontsize',12);
uitext(u,lines1,lines2,lines3,lines4,lines5,lines6);


if(numGoodMotifs > 0)
    %Display Example Traces:
    displayCellMotif(cell, cell.nMotifBest, numVert, numHorz, 5, 7, false);
    displayCellMotif(cell, cell.nMotifWorst, numVert, numHorz, 6, 8, false);
    displayCellMotif(cell, cell.nMotifBest2, numVert, numHorz, 9, 11, false);
    displayCellMotif(cell, cell.nMotifWorst2, numVert, numHorz, 10, 12, false);  
 
    %Display a Raster
    
    subplot(numVert,numHorz,[13:16]);
    cellRaster(cell);

    blur = .005;
    [R, commonStartTime, commonEndTime, motifMap] = cellMotifPairWiseCorrelations(cell, false, false, blur);
    
    %Pair wise correlation histogram
    subplot(numVert,numHorz,[18,20]);
    hist(reshape(R,1,prod(size(R))),[-1:.01:1]);
    mu = mean(R(find(R<1)));
    h = line([mu,mu],ylim); 
    set(h,'Color','red');
    std = std(R(find(R<1)));
    h = line([mu-std,mu+std],[5,5]);set(h,'Color','red'); 
    xlim([-1,1]);
    legend(['mu = ', num2str(mu)]); 
    
    %Pair wise correlation matrix.
    subplot(numVert,numHorz,[17,19]);
    imagesc(R);
    xlabel(['Pairwise Corr. Gauss=',num2str(blur),' startT=',num2str(commonStartTime),' endT=',num2str(commonEndTime)]);
    axis square;
end
    
    
    
    
    
    
     
