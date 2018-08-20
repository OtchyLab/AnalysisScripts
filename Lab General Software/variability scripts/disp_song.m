function [song,filename]=disp_song(n, all_files, file_order)

    filename=all_files(file_order(n-2)).name
    song=wavread(filename);
    specgram1(song,512,44100,400,350);ylim([0 10000]);
end