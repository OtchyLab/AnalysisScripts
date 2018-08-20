function multipleCellRaster(cellNums, cellDir, startTime, endTime)

colors = {'black','red','blue','green','magenta'};
numRows = 0;
nCell = 0;
for(cellNum = cellNums)
    nCell = nCell + 1;
    cell = loadCell(cellNum,cellDir)
    color = colors{mod(nCell,length(colors))+1};
    numRows = cellRasterWithLines(cell, [], 1, length(cell.motifData), [], startTime, endTime, color, numRows);
end

