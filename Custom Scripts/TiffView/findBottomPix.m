function [row,col] = findBottomPix(binImage)
    % bin image is a 2d array of doubles with binary coloring
    [numRows, numCols] = size(binImage);
    found = false;
    curRow = numRows; 
    while found == false
        for i = 1:numCols
            if (binImage(curRow, i) == 0)
                row = curRow;
                col = i;
                found = true;
                break;
            end
        end
        curRow = curRow-1;
    end
end