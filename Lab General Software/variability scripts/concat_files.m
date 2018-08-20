function concat_files
Fs=44100;
[filename,dirpath]=uigetfile('*.*');
old_dir=cd;
cd (dirpath);
all_files=dir(dirpath);
[file_order,file_no]=sortfiles(all_files, filename);
filename=all_files(file_order(1)).name;
full_song=0;
buffer=zeros(2000,1);
for i=3:length(all_files)
    filename=all_files(file_order(i-2)).name;
    song=wavread(filename);
    norm_factor=std(song);
    song=song*(1/(20*norm_factor));
    full_song=[full_song;buffer;song];
end
[FileName,PathName] = uiputfile;
cd (PathName);    
wavwrite(full_song,Fs,FileName);
cd (old_dir);
