for i = 1:length(AllCells)
    cell = loadCell(AllCells(i));
    [R, startTime, endTime, motifMap] = cellMotifPairWiseCorrelations(cell, [], false, true, .02)

% for i = Included
%     cell = loadCell(i);
%     rawCount(i) = length(cell.motifData);
%     goodCount(i) = 0;
%     for(nMotif = 1:rawCount(i))
%         if(isfield(cell.motifData{nMotif},'spikeNdx'))
%             goodCount(i) = goodCount(i) +1;
%         end
%     end
% end


% 
% h1 = figure; 
% edges = 10.^([-4:.04:1]);
% n1 = histc(isi_good_song,edges);
% n1 = n1/sum(isi_good_song);
% h = bar(edges(1:end), n1);
% xlabel('ISI in seconds');
% ylabel('Count');
% title(['Normalized Song']);
% h1  = get(h1,'Children');
% set(h,'LineWidth'
% set(h1,'XScale','log');
% 
% h2 = figure; 
% n2 = histc(isi_good_nonsong,edges);
% n2 = n2/sum(isi_good_nonsong);
% bar(edges(1:end), n2);
% %bar(edges, n2,'histc');
% xlabel('ISI in seconds');
% ylabel('Count');
% title(['Normalized Non Song']);
% h2  = get(h2,'Children');
% set(h2,'XScale','log');




% isi_good_song = [];
% isi_good_nonsong = [];
% for(i = [1,2,4,5,7])
%     isi_good_song = [isi_good_song, isiSong{i}];
%     isi_good_nonsong = [isi_good_nonsong, isiNonSong{i}];
% end
% 
% length(isi_good_song)
% length(isi_good_nonsong)
% 
% 
% figure; 
% hist(isi_good_song,0:.0001:.1);
% xlabel('ISI in seconds');
% ylabel('Count');
% title(['good']);
% 
% figure; 
% hist(isi_good_nonsong,0:.0001:.1);
% xlabel('ISI in seconds');
% ylabel('Count');
% title(['good non']);

goodsong = length(find(isi_good_song<.005)) / length(isi_good_song)
goodnonsong = length(find(isi_good_nonsong<.005)) / length(isi_good_nonsong)

% all_isiSong = []
% all_isiNonSong = []
% for(nCell = 1:length(Collision))
%     [isiSong{nCell}, isiNonSong{nCell}] = compareISISingingToNonSinging(loadCell(Collision(nCell)));
%     all_isiSong = [all_isiSong, isiSong{nCell}];
%     all_isiNonSong = [all_isiNonSong, isiNonSong{nCell}];
% end

% for(nCell = 1:length(Collision))
%     fracBurstSong(nCell) = length(find(isiSong{nCell}<.010)) / length(isiSong{nCell});
%     fracBurstNonSong(nCell) = length(find(isiNonSong{nCell}<.010)) / length(isiNonSong{nCell});
%     firingRateSong(nCell) = length(isiSong{nCell}) / sum(isiSong{nCell});
%     firingRateNonSong(nCell) = length(isiNonSong{nCell}) / sum(isiNonSong{nCell});
% end
% 
% allFracBurstSong = length(find(all_isiSong<.010)) / length(all_isiSong)
% allFracBurstNonSong = length(find(all_isiNonSong<.010)) / length(all_isiNonSong)
% mean(fracBurstSong)
% mean(fracBurstNonSong)

 