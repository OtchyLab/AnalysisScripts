function CellFileConvert
%This function is to convert the old, BDT-file based cell standard into
%the new, DAT-file based cell standard. The main effect of this script is
%to cross reference the filenum of the old cell with newer,
%timestamped ID method. Some of the information contained in the old
%cell file is reformatted or discarded; additional info is added.

%Pre-allocate main variables
old_presentCell = [];
new_presentCell = [];
DATnames = [];
mother = 'V:';

%Select annotations to add compile and combine
[oldFiles,oldPathName] = uigetfile([mother, '\Single Units\*.mat'],'Select cell file to convert.');

% %Get the location of the original DAT-files to cross-reference
% DATpath = uigetdir('V:\Pur188\','Select location of the DAT files.');
% 
% %Get the location to save the output file
% [newFile,newPathName] = uiputfile(['V:\Pur188\Cells\' oldFiles],'Where should the new cell file be saved?');

%parse the cell name for suggestion
s = regexp(oldFiles, '_', 'split');
birdName = s{1};
dates = s{2};
depth = s{3};

d = dates(end-1:end);
m = dates(end-3:end-2);
y = dates(1:end-4);

%Convert to 4-digit year
if numel(y) == 2
    y = ['20', y];
end

%Set location of the original DATs
DATpath = [mother filesep birdName filesep y '-' m '-' d filesep];

%Set output location
newPathName = [mother filesep birdName filesep 'Cells' filesep];
newFile = [birdName '_' y m d '_' depth '_cell.mat'];

if ~exist(newPathName)
    mkdir(newPathName)
end

%Load cell file and stack all in a common variable
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
    
    if exist('presentCell')
        
        old_presentCell = [old_presentCell, presentCell];
    else
        disp(['Unexpected content in file ' oldFiles{i} '. Check format and try again.'])
        return
    end
    clear('presentCell');
end

%Index all of the DAT files in the chosen folder
% temp = dir([DATpath '\*.DAT']);
DATfiles = dir([DATpath '\*.dat']);
% DATfiles = [DATfiles; temp];
for i = 1:length(DATfiles)
    DATnames{i} = DATfiles(i).name(19:end-4);
end

%Loop through the old annotation file by file
index = 1;
for i = 1:length(old_presentCell)
    %Parse the BDT-file name for the date and time info
    %temp = regexp(old_keys{i},'_','split');
    if ~isempty(old_presentCell(i).filename) && ~strcmp(old_presentCell(i).filename, '0')
        temp = regexp(old_presentCell(i).filename,'_','split');
        Y = temp{end}(1:4); M = temp{end}(5:6); D = temp{end}(7:8);
        h = temp{end}(10:11); m = temp{end}(12:13); s = temp{end}(14:15);
    
        %Create the matching string template
        template = [Y '_' M '_' D '_' h '_' m '_' s];
    
        %Search for the matching string
        pntr = strmatch(template,DATnames,'exact');
        if isempty(pntr)
            %Report no match found
            disp(['No match found for ' old_presentCell(i).filename])
        elseif length(pntr)>1
            %Report more than one match found
            disp(['More than one entry found for ' old_presentCell(i).filename])
        else
            %Update the record in the new cell variable
             record.birdname = old_presentCell(i).birdname;
             record.channel = old_presentCell(i).channel;
             
             record.cell_no = old_presentCell(i).cell_no;

             temp = regexp(DATfiles(pntr).name,'_','split');
             record.filenum = str2num(temp{2});

%              if strcmp(old_presentCell(i).filename(1:5),'times')
%                 record.filename = ['times_' DATfiles(pntr).name];
%              else
                 record.filename = DATfiles(pntr).name;
%              end

             %Formatting check
             if ~ischar(record.channel)
                 record.channel = num2str(record.channel);
             end

             %Empirically determined on 2/19/16; use this code to convert cell files that have both the 
             % 44.1/44.15 kHz difference and the MUX offsets.
             c = 0.0011; % based on 50 samples @ 44150Hz             
             if strcmp(record.channel,'1')
                 offset = 1*c;
             elseif strcmp(record.channel,'2')
                 offset = 2*c;
             elseif strcmp(record.channel,'3')
                 offset = 3*c;
             elseif strcmp(record.channel,'4')
                 offset = 4*c;
             end
            scale = 0.998867497168743;
            record.spikes = scale * [old_presentCell(i).spikes(old_presentCell(i).spikes < 1) - offset; old_presentCell(i).spikes(old_presentCell(i).spikes >= 1)];
            record.noise = scale * [old_presentCell(i).noise(old_presentCell(i).noise < 1) - offset; old_presentCell(i).noise(old_presentCell(i).noise >= 1)];
            record.spont = scale * [old_presentCell(i).spont(old_presentCell(i).spont < 1) - offset; old_presentCell(i).spont(old_presentCell(i).spont >= 1)];


            %Empirically determined on 2/19/16; seems principally due to the 44.1/44.15 kHz difference
            %Looks to work well for all channels if the only issue is trhe sampling rate miscalculation
%              offset = 0; %Null out for now...
%              scale = 0.998867497168743; 
% 
%              record.spikes = (scale*old_presentCell(i).spikes) - offset;
%              record.noise = (scale*old_presentCell(i).noise) - offset;
%              record.spont = (scale*old_presentCell(i).spont) - offset;

            %Set constant sampling rate
            record.fs = 44150;
            
            %Copy to the new structures
            new_presentCell = [new_presentCell; record];

            %Update output index
            index = index + 1;
        end
    end
end
presentCell = new_presentCell';

%Save the output file
save([newPathName,newFile],'presentCell', '-v7');
disp(['New cell file created at ' newPathName newFile]);










