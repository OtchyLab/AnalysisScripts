folder = uigetdir('C:\User\Tim\Desktop');

cd(folder)
files = dir('*.mat');;
%files.name = 'Sil167_1835cell_112dph.mat';
numFiles = length(files);

for i=1:numFiles
    curFile = char(files(i).name);
    load(curFile,'audio','spikes');
    t = curFile(1:end-4);
    numRend = length(audio);
    mkdir(t);
    for j=1:numRend
       nameWAV = [t '_r' num2str(j) '.wav'];
       nameTXT = [t '_r' num2str(j) '.txt'];
       wavAud = audio{j}./max(audio{j});
       wavwrite(wavAud,44150,[folder '\' t '\' nameWAV]);
       spiketxt = spikes{j};
       %save([folder '\' t '\' nameTXT], 'spiketxt','-ascii', '-double');
       fid = fopen([folder '\' t '\' nameTXT],'wt');  % Note the 'wt' for writing in text mode
       fprintf(fid,'%f\n',spiketxt);  % The format string is applied to each element of a
       fclose(fid);
    end
end
