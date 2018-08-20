function annotationResampler
%This function is to adjust the annotation syllable times to account for the mis-sampling of the raw data in the pre-2011
%era. The code is easily adapted to account for arbitrary up/down sampling, but is now writting for the 44100 ---> 44150Hz
%transition.

%Sampling conversion
origFS = 44100; %in Hz
newFS = 44150;

%Conversion ratio
r = origFS/newFS;

%Ask the user to choose the annotation file to convert
[oldFiles,oldPathName] = uigetfile('V:\*.mat','Select annotation file to resample.','MultiSelect','on');
%oldFile = 'Sil167_20100216_1870u_annotation.mat';
%oldPathName = 'V:\Single Units\Sil 167 Sorted Cells\Sil167 021610 111dph\RA record\';

%Get the location to save the output file
%[newFile,newPathName] = uiputfile('V:\*.mat','Where should the new annotation be saved?');

%Cycle through the selected annotation files
for i = 1:numel(oldFiles)
    
    %Load the annotaion
    keys = []; elements =[];
    load([oldPathName oldFiles{i}]);

    %Check if properly loaded
    if ~iscell(keys) || ~iscell(elements)
        display(['The file load is fucked up... check file #' num2str(i)])
        beep; return
    end

    %Rescale time
    for j = 1:numel(elements)
        %Retain the old values
        old = elements{j}.segFileStartTimes;

        %Scale the relative starts and ends by the FS ratio
        newS = r.*elements{j}.segFileStartTimes;
        newE = r.*elements{j}.segFileEndTimes;

        %Temporal difference between old and new
        delta = newS-old;

        %Adjust absolute times
        absS = elements{j}.segAbsStartTimes + delta;

        %Copy back the values to the elements structure
        elements{j}.segFileStartTimes = newS;
        elements{j}.segFileEndTimes = newE;
        elements{j}.segAbsStartTimes = absS;

    end

    %Create the new filename
    base = oldFiles{i}(1:end-4);
    newFile = [base, 'RS.mat'];
    
    %Save keys and elements to the specified location
    save([oldPathName newFile], 'keys', 'elements', '-v7');
    

end

