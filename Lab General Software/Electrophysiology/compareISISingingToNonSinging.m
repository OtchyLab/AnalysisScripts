function [isiSong, isiNonSong] = compareISISingingToNonSinging(cell)

h = figure;
subplot(2,1,1);
isiSong = getCellISIDuringSong(cell);
hist(isiSong,0:.0001:.1);
xlabel('ISI in seconds');
ylabel('Count');
title(['ISI distribution during Song Cell#',num2str(cell.num)]);
    
subplot(2,1,2);
isiNonSong = getCellISINotSinging(cell);
hist(isiNonSong,0:.0001:.1);
xlabel('ISI in seconds');
ylabel('Count');
title(['ISI distribution when non singing Cell#',num2str(cell.num)]);