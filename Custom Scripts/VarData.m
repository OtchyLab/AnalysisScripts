function VarData(end_date, earliest_date, interval)
%Top level function call for batch stripping the audio data into single
%syllable files and and saving in an orderly way.  Based on the SylSim
%routine.

%Ask user which file to find the dated folders of source data
%parentsource = uigetdir(pwd,'Source Folder');
parentsource = 'G:\Manually Backed Up Audio\Rd278';

%Ask user where to store syllables
%parentdest = uigetdir(pwd,'Destination Folder');
parentdest = 'C:\Users\Tim\Desktop\Rd278 Chopped Syls';

date = end_date;

while datenum(date, 'yyyy-mm-dd') >= datenum(earliest_date, 'yyyy-mm-dd')
    songsource = [parentsource '\' date];
    while exist(songsource)~=7
        date = datestr(addtodate(datenum(date), -1, 'day'), 'yyyy-mm-dd');
        songsource = [parentsource '\' date];
    end
    if isdir(songsource)
        dest = [parentdest '\' date];
        if ~isdir(dest)
            mkdir(dest);
        end
        SylParam = []; DataSet = []; PCA_set = []; sim = [];
        SylParam = genvarname(['SylParam_', datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd')]);
        DataSet = genvarname(['DataSet_', datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd')]);
        RawDataSet = genvarname(['RawDataSet_', datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd')]);
        eval(['[' SylParam ', ' DataSet ', ' RawDataSet ']=SylSim(songsource, dest);']);
        cd(dest);
        eval(['save ' date ' ' SylParam ' ' DataSet ' ' RawDataSet ';']);
    else
    end
    date = datestr(addtodate(datenum(date), -1*interval, 'day'), 'yyyy-mm-dd');
end
end
