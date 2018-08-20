%script M-file "readwavefile.m". This reads in a wavefile in a specified
%directory

[filename,dirpath]=uigetfile('*.*');
old_dir=cd;
cd (dirpath);
all_files=dir(dirpath);
[file_order,file_no]=sortfiles(all_files, filename);
n=file_no;
%while (n<=length(all_files))%&& (stop button has not been pressed)
    filename=all_files(file_order(n-2)).name
    song=wavread(filename);
    display_spect;
    axes(a);
    specgram1(song,512,44100,400,350);
 
    n=n+1;
%end
cd(old_dir);
