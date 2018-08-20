function AnnotationConvert
%This function is to convert the old, BDT-file based annotation standard into
%the new, DAT-file based annotation standard. The main effect of this script is
%to cross reference the filenum of the old annotation with newer,
%timestamped ID method. Some of the information contained in the old
%annotation file is reformatted or discarded; additional info is added.

%Pre-allocate main variables
old_elements = [];
new_elements = [];
old_keys = [];
new_keys = [];
DATnames = [];

%Select annotations to add compile and combine
[oldFiles,oldPathName] = uigetfile('V:\Single Units\Sil 167 Sorted Cells (WT)\*.mat','Select all annotation files to include.','MultiSelect','on');
%oldFiles = 'Sil167_20100216_1870u_annotation.mat';
%oldPathName = 'V:\Single Units\Sil 167 Sorted Cells\Sil167 021610 111dph\RA record\';

%Get the location of the original DAT-files to cross-reference
DATpath = uigetdir('V:\Sil167\','Select location of the DAT files.');
%DATpath = 'V:\Sil167\2010-02-16';

%Get the location to save the output file
[newFile,newPathName] = uiputfile('V:\Sil167\*.mat','Where should the new annotation be saved?');
%newFile = 'testAnnot.mat';
%newPathName = 'V:\Pur626\';

%Load annotations and stack all in a common variable
if iscell(oldFiles)
    numAnnots = size(oldFiles,2);
else
    numAnnots = size(oldFiles,1);
end
for i = 1:numAnnots
    if iscell(oldFiles)
        load([oldPathName oldFiles{i}])
    else
        load([oldPathName oldFiles])
    end
    
    if exist('elements') && exist('keys')
        old_elements = [old_elements, elements];
        old_keys = [old_keys, keys];
    else
        disp(['Unexpected content in file ' oldFiles{i} '. Check format and try again.'])
        return
    end
    clear('elements','keys');
end

%Index all of the DAT files in the chosen folder
temp = dir([DATpath '*.DAT']);
DATfiles = dir([DATpath '\*.dat']);
DATfiles = [DATfiles; temp];
for i = 1:length(DATfiles)
    DATnames{i} = DATfiles(i).name(19:end-4);
end

%Loop through the old annotation file by file
index = 1;
for i = 1:length(old_keys)
    %Parse the BDT-file name for the date and time info
    temp = regexp(old_keys{i},'_','split');
    Y = temp{3}(1:4); M = temp{3}(5:6); D = temp{3}(7:8);
    h = temp{3}(10:11); m = temp{3}(12:13); s = temp{3}(14:15);
    
    %Create the matching string template
    template = [Y '_' M '_' D '_' h '_' m '_' s];
    
    %Search for the matching string
    pntr = strmatch(template,DATnames,'exact');
    if isempty(pntr)
        %Report no match found
        disp(['No match found for ' old_keys(i)])
    elseif length(pntr)>1
        %Report more than one match found
        disp(['More than one entry found for ' old_keys(i)])
    else
        %Update the record in the new annotation variable
        record.exper.birdname = old_elements{i}.exper.birdname;
        record.exper.expername = [old_elements{i}.exper.birddesc ' ' old_elements{i}.exper.experdesc];
        record.exper.desiredInSampRate = 44150;
        record.exper.audioCh = old_elements{i}.exper.audioCh;
        record.exper.sigCh = old_elements{i}.exper.sigCh;
        record.exper.datecreated = datestr(now,'dd-mmm-yyyy');
        record.exper.researcher = 'TMO';
        
        temp = regexp(DATfiles(pntr).name,'_','split');
        record.filenum = str2num(temp{2});
        record.segAbsStartTimes = old_elements{i}.segAbsStartTimes;
        record.segFileStartTimes = old_elements{i}.segFileStartTimes;
        record.segFileEndTimes = old_elements{i}.segFileEndTimes;
        record.segType = old_elements{i}.segType;
        record.fs = 44150;
        record.drugstatus = old_elements{i}.drugstatus;
        record.directstatus = 'Undirected';
        
        %Copy to the new structures
        new_elements{index} = record;
        new_keys{index} = DATfiles(pntr).name;
        
        %Update output index
        index = index + 1;
    end
end
[keys,IX] = sort(new_keys);
elements = new_elements(IX);

%Save the output file
save([newPathName,newFile],'elements','keys');
disp(['New annotation file created at ' newPathName newFile]);










