function viewCellsISI(cellList)

h = figure;
numPerFigure = 2;
%Show numPerFigure cells at a time in subplots 1 through 5.
for(nCell = 1:length(cellList))
    subplot(numPerFigure,1,mod(nCell-1,numPerFigure) + 1);
    cell = loadCell(cellList(nCell));
    isi = getCellISIDuringSong(cell);
    hist(isi,0:.0001:.1);
    xlabel('ISI in seconds');
    ylabel('Count');
    title(['ISI distribution Cell#',num2str(cell.num)]);
    
    if(mod(nCell,numPerFigure) == 0)
        figure(h);
        set(h,'CurrentCharacter','z');
        while(true)
            pause(.05);
            char = get(h,'CurrentCharacter');
            if(char == ' ')
                break;
            end
        end
    end    
end