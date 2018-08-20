function [R, Rvalues, allRvalues, muR, stdR, countR] = correlationsForFigure(cellList, gaussWidth, bRotate)

nCell = 1;
allRvalues = []
for cellNum = cellList
    cell = loadCell(cellNum)
    R{nCell} = cellMotifPairWiseCorrelations(cell, [], false, bRotate, gaussWidth);
    Rvalues{nCell} = R{nCell}(find(R{nCell}<1));

    if(bRotate)
        for(nEnsemble = 2:100)
            R{nCell} = cellMotifPairWiseCorrelations(cell, [], false, bRotate, gaussWidth);
            Rvalues{nCell} = [Rvalues{nCell}; R{nCell}(find(R{nCell}<1))];          
        end
    end
    
    allRvalues = [allRvalues; Rvalues{nCell}];
    muR(nCell) = mean(Rvalues{nCell});
    stdR(nCell) = std(Rvalues{nCell});
    countR(nCell) = length(Rvalues{nCell});
    disp(nCell);
    nCell = nCell + 1;
end

muAll = mean(allRvalues);
stdAll = std(allRvalues);
muCell = mean(muR);
