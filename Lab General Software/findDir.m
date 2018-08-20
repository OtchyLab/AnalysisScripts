function findDir

[dirpath]=uigetdir;
cd(dirpath);
directory=dir;
if isempty(directory)
    warndlg('no directories found')
    return;
end
for i=1:length (directory)
    pathname = [dirpath,'\', directory(i).name];
    filteraudio(pathname);
end
