function cell = chooseBestAndWorstTrace(cell, startMotif)

cell = loadCell(cell.num);

if(~exist('startMotif'))
    startMotif = 1;
end

opt = '~';
nSecondWorst = 0;
nSecondBest = 0;
nWorst = 0;
nBest = 0;
h = figure(1006);
for(nMotif = startMotif:length(cell.motifData))
    motifData = cell.motifData{nMotif};
    if(isfield(motifData,'spikeNdx'))
        displayCellMotif(cell,nMotif,4,1,1,2);
        title(num2str(nMotif));
        figure(h);
        while(true)
            pause(.05);
            opt = get(h,'CurrentCharacter');
            figure(h);
            set(h,'CurrentCharacter','~')
            if(opt == '.')
                break;
                figure(h);
            elseif(opt == 'q')
                break;               
            elseif(opt == 'b')
                nSecondBest = nBest;
                nBest = nMotif;
                subplot(4,1,3);
                displayMotifTrace(cell,nBest);
                title(['Cell ', num2str(cell.num), ' Motif ',num2str(nBest)]);
                figure(h);
            elseif(opt == 'w')
                nSecondWorst = nWorst;
                nWorst = nMotif;
                subplot(4,1,4);
                displayMotifTrace(cell,nWorst);
                title(['Cell ', num2str(cell.num), ' Motif ',num2str(nWorst)]);
                figure(h);
            elseif(opt == 'c')
                if(~isfield(cell,'comments'))
                    nComment = 1;
                else
                    nComment = length(cell.comments) + 1;
                end
                cell.comments{nComment} = input('Comment: ','s');  
                figure(h);
            end
        end
    end
    if(opt == 'q')
        break;
    end
end

cell.nMotifBest = nBest;
cell.nMotifBest2 = nSecondBest;
cell.nMotifWorst = nWorst;
cell.nMotifWorst2 = nSecondWorst;
saveCell(cell);
    