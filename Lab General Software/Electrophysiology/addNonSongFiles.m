function addNonSongFiles(cellList)

for(nCell = 1:length(cellList))
    cell = loadCell(cellList(nCell));
    createNonSongCellFiles(cell);
end