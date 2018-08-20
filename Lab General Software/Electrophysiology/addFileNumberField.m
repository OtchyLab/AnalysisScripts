function addFileNumberField(cellList)

for(nCell = 1:length(cellList))
    cell = loadCell(cellList(nCell));
    cell.fileNumbers = input(['Enter file numbers labnotebook cell #',num2str(cell.labNotebookNum)]);
    saveCell(cell);
end