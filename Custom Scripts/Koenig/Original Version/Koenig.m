function varargout = Koenig(varargin)
% KOENIG MATLAB code for Koenig.fig
%      KOENIG, by itself, creates a new KOENIG or raises the existing
%      singleton*.
%
%      H = KOENIG returns the handle to a new KOENIG or the handle to
%      the existing singleton*.
%
%      KOENIG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KOENIG.M with the given input arguments.
%
%      KOENIG('Property','Value',...) creates a new KOENIG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Koenig_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Koenig_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Koenig

% Last Modified by GUIDE v2.5 14-Nov-2012 17:21:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Koenig_OpeningFcn, ...
                   'gui_OutputFcn',  @Koenig_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

function Koenig_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Koenig (see VARARGIN)

% Choose default command line output for Koenig
handles.output = hObject;

handles.annotation = mhashtable;

%set(handles.popupSequence,'Value',1);

handles.fs=44150;
handles.dirname =[];
handles.annotationFilename =[];
handles.elements = [];
handles.totFilelist = [];
handles.filelist =[];
handles.dataset =[];

%Clear axes
set(handles.axes_timeDistr,'XTick',[],'YTick',[]);
set(handles.axes_sylDistr,'XTick',[],'YTick',[]);
set(handles.axes_transMatrix,'XTick',[],'YTick',[]);
set(handles.axes5,'XTick',[],'YTick',[]);
set(handles.axes6,'XTick',[],'YTick',[]);
set(handles.axes_SylGapDistr,'XTick',[],'YTick',[]);
set(handles.axes_SylIntDur,'XTick',[],'YTick',[]);
set(handles.axes13,'XTick',[],'YTick',[]);
set(handles.axes_rhythm,'XTick',[],'YTick',[]);
set(handles.axes_featStats,'XTick',[],'YTick',[]);
set(handles.axes_featVfeat,'XTick',[],'YTick',[]);
set(handles.axes_simStats,'XTick',[],'YTick',[]);

%Set initial values for feature popup boxes
handles.featureSet = {'Duration';'Amp Mod';'Freq Mod';'Entropy';'1/2 Amp';'Grav Cent';...
    'Pitch Good'; 'Pitch'};
set(handles.popup_featStats,'String',handles.featureSet);
set(handles.popup_feat1,'String',handles.featureSet);
set(handles.popup_feat2,'String',handles.featureSet);

% Update handles structure
guidata(hObject, handles);

function varargout = Koenig_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Independent, in-GUI functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles] = clearAllDataAxes(handles)
% Clears out all of the axes that show changeable data
axes(handles.axes_timeDistr)
cla;
axes(handles.axes_sylDistr)
cla;
axes(handles.axes_transMatrix)
cla;
axes(handles.axes5)
cla;
axes(handles.axes6)
cla;
axes(handles.axes_SylGapDistr)
cla;
axes(handles.axes_SylIntDur)
cla;
axes(handles.axes13)
cla;
axes(handles.axes_rhythm)
cla;
axes(handles.axes_featStats)
cla;
axes(handles.axes_featVfeat)
cla;
axes(handles.axes_simStats)
cla;

function [xRefFiles] = XrefAnnotData(totalFiles, keys)
% This function takes in a the sorted string list "keys" which contains all
% of the .dat filenames that are required by the annotation files and
% "totalFiles" which is a structure list of all the .dat files contained
% in the loaded data folders.  "xRefFiles" will return the subset of
% totalFiles that correspond to each row of "keys".  In the event that the
% matching file cannot be found, load a dummy record in it's place

%Preallocate
fileInd = zeros(length(keys),1);
flnm = {};

%Define Dummy Record for gaps
dummy.name = 'DUMMY';
dummy.date = 'DUMMY';
dummy.bytes = 6969;
dummy.isdir = false;
dummy.datenum = 6969;
dummy.dirname = 'DUMMY';

%Extract out the filenames of all .dat files
for i = 1:length(totalFiles)
    flnm{i} = totalFiles(i).name;   
end

for i = 1:length(keys)
    index = strmatch(keys(i),flnm);
    if length(index)>1
        xRefFiles(i) = totalFiles(index(1));
    elseif isempty(index)
        xRefFiles(i) = dummy;
    else
        xRefFiles(i) = totalFiles(index);
    end
end

function [filter_Ind] = filterRecords(handles)
%The purpose of this function is to take in the complete handles.filelist
%(which contains all the files that are in both the data folders and in the
%annotations) and return a new list that has been appropriately parsed by
%the selections made by the user.

%Create index for drug status
if strcmp(handles.filter_DrugStat,'All')
    drugInd = ones(length(handles.drugstatus),1)';
else
    drugInd = strcmp(handles.filter_DrugStat,handles.drugstatus);
end

%Create index for directed status
if strcmp(handles.filter_DirectStat,'All')
    directInd = ones(length(handles.directstatus),1)';
else
    directInd = strcmp(handles.filter_DirectStat,handles.directstatus);
end

%Create index for the file number limits
recInd = (handles.filenums>=handles.filter_fromRec & handles.filenums<=handles.filter_toRec);

%Combine all indices to master filter index
filter_Ind = bitand(bitand(drugInd,directInd),recInd);

function [allSeqNumbers startSeq]=sortSequences(handles,numberSyll)
%make a cell array with all the song sequences

[allSeq startSeq counter]=getSequences(handles,numberSyll);
if (allSeq==0)
    warndlg('There are no sequences in the specified range.');
    uiwait;
    return;
end
%go through this cell array and build up an array of sequences
%numberSyll long
for i=1:numberSyll
    multiplier(i)=10^(numberSyll-i);
end
allSeqNumbers=[];
for i=1:counter-1
    allSeqNumbers(i)=dot(allSeq(i,:),multiplier);
end

function [allSeq startSeq counter]=getSequences(handles,numberSyll)
%Gets all the sequences containing numSyll.  Note that this is entirely
%dependent on the current filter settings.

%Reset placeholders
counter=1;
song_interruptions=[];

%Restrict sequence search to those files filtered by the user
filt_elements = handles.elements(handles.filtInd == 1);

%Find all pauses greater than 100ms to determine boundaries of song motifs
for i = 1:length(filt_elements)
    for j = 1:length(filt_elements{i}.segFileStartTimes)-1
        pauses(j) = filt_elements{i}.segFileStartTimes(j+1)-filt_elements{i}.segFileEndTimes(j);
        song_interruptions=find(pauses>0.1);
    end

    if length(filt_elements{i}.segFileStartTimes)>=numberSyll
        for j = 1:length(filt_elements{i}.segFileStartTimes)-numberSyll+1
            current_seq(1:numberSyll)=filt_elements{i}.segType(j:j+numberSyll-1);
            if (isempty(find(current_seq<1 | current_seq>9)) && isempty(find(song_interruptions>j-1 & song_interruptions<j+numberSyll-1)))
                allSeq(counter,1:numberSyll)=current_seq;
                startSeq(counter,:)=[i,j];%index (index in annotation file,number of starting syllable of sequences)
                counter=counter+1;
            end
        end
    end
end

if ~exist('allSeq')
    allSeq=0;
    startSeq=0;
    counter=0;
end

function [syll motif audio] = getSubSet(handles,buffer)
%Given all of the previous filtering and selecting, this function will
%parse the data from the dataset, perform the requested alignments, and
%display it on both the audio and neuro axes.

%Parse data from the dataset into two matrices (audio and neuro)
%Filter the annotation and dataset
filt_elements = handles.elements(handles.filtInd==1);
filt_dataset = handles.dataset(handles.filtInd==1);

for i=1:size(handles.chosenStartSeq,1)
    %Grab the start and end times of each syllable in the motif 
    syll(1:2:handles.numSylls*2-1,i)=filt_elements{handles.chosenStartSeq(i,1)}.segFileStartTimes(handles.chosenStartSeq(i,2):handles.chosenStartSeq(i,2)+handles.numSylls-1);
    syll(2:2:handles.numSylls*2,i)=filt_elements{handles.chosenStartSeq(i,1)}.segFileEndTimes(handles.chosenStartSeq(i,2):handles.chosenStartSeq(i,2)+handles.numSylls-1);
    
    %Grab the start and end times of each motif
    motif(1,i)= filt_elements{handles.chosenStartSeq(i,1)}.segFileStartTimes(handles.chosenStartSeq(i,2));
    motif(2,i)= filt_elements{handles.chosenStartSeq(i,1)}.segFileEndTimes(handles.chosenStartSeq(i,2)+handles.numSylls-1);

    %Calculate the boundaries for snipping
    startT = max(0,floor(((motif(1,i)-buffer)*handles.fs)));
    endT = ceil(((motif(2,i)+buffer)*handles.fs));
    
    %Grab the segment of data that corresponds to each selected motif
    audio{i} = filt_dataset{handles.chosenStartSeq(i,1)}(1,startT:endT);
end

function [subAudio subElements] = getWorkset(handles)
if get(handles.check_useAllSyls,'Value')
    subAudio = handles.dataset(handles.filtInd==1);
    subElements = handles.elements(handles.filtInd==1);
else
    %This should be where the motif-restricted workset is created... not
    %obvious to me how this should be done yet.
end

function [transMatrix,sylTypes] = SequenceMatrixSE(elements)
%go through the annotations and whenever there is a stretch of contiguous
%sequence assign the transition to the matrix. display the figure.

%Reset placeholders
counter=1;
song_interruptions=[];
numberSyll=2;

%Find all pauses greater than 100ms to determine boundaries of song motifs
for i = 1:length(elements)
    pauses = [];
    for j = 1:length(elements{i}.segFileStartTimes)-1
        pauses(j) = elements{i}.segFileStartTimes(j+1)-elements{i}.segFileEndTimes(j);
        song_interruptions=find(pauses>0.1);
    end

    if length(elements{i}.segFileStartTimes)>=numberSyll
        %Reshape the syllable sequence for the whole file to reflect the
        %start and ends of song bouts
        segTypeExt = elements{i}.segType;
        boutbreak = [99;-99];
        for j = 1:length(song_interruptions)
            ind = song_interruptions(j)+(2*(j-1));
            pre = segTypeExt(1:ind);
            post = segTypeExt((ind+1):end);
            segTypeExt = [pre; boutbreak; post];
        end
        segTypeExt = [-99; segTypeExt; 99];
        
        for j = 1:length(segTypeExt)-numberSyll+1
            current_seq(1:numberSyll) = segTypeExt(j:j+numberSyll-1);
            %if (isempty(find(current_seq<1 | current_seq>9)) && isempty(find(song_interruptions>j-1 & song_interruptions<j+numberSyll-1)))
            if (isempty(find(current_seq==-1 | current_seq>99)) && ~isequal(current_seq,boutbreak'))
                allSeq(counter,1:numberSyll)=current_seq;
                startSeq(counter,:)=[i,j];
                counter=counter+1;
            end
        end
    end
end

%numberSyllTypes=max(max(allSeq));
numberSyllTypes=length(unique(allSeq));
sylTypes = sort(unique(allSeq));
transMatrix=zeros(numberSyllTypes,numberSyllTypes);   
for i=1:counter-1
    from=find(sylTypes==allSeq(i,1));
    to=find(sylTypes==allSeq(i,2));
    transMatrix(from,to)=transMatrix(from,to)+1;
end

function [transMatrix,sylTypes] = SequenceMatrix(elements)
%go through the annotations and whenever there is a stretch of contiguous
%sequence assign the transition to the matrix. display the figure.

%Reset placeholders
counter=1;
song_interruptions=[];
numberSyll=2;

%Find all pauses greater than 100ms to determine boundaries of song motifs
for i = 1:length(elements)
    for j = 1:length(elements{i}.segFileStartTimes)-1
        pauses(j) = elements{i}.segFileStartTimes(j+1)-elements{i}.segFileEndTimes(j);
        song_interruptions=find(pauses>0.1);
    end

    if length(elements{i}.segFileStartTimes)>=numberSyll
        for j = 1:length(elements{i}.segFileStartTimes)-numberSyll+1
            current_seq(1:numberSyll)=elements{i}.segType(j:j+numberSyll-1);
            if (isempty(find(current_seq<1 | current_seq>9)) && isempty(find(song_interruptions>j-1 & song_interruptions<j+numberSyll-1)))
                allSeq(counter,1:numberSyll)=current_seq;
                startSeq(counter,:)=[i,j];%index (index in annotation file,number of starting syllable of sequences)
                counter=counter+1;
            end
        end
    end
end

numberSyllTypes=max(max(allSeq));
transMatrix=zeros(numberSyllTypes,numberSyllTypes);   
for i=1:counter-1
    from=allSeq(i,1);
    to=allSeq(i,2);
    transMatrix(from,to)=transMatrix(from,to)+1;
end

sylTypes = unique(allSeq);

function [sylDurs,gapDurs,sylDursbyType,sylTypes] = SylGapDistr(elements)
%Initialize vars
sylLabels = [];
sylDurs = [];
gapLabels = [];
gapDurs = [];
starts = [];
ends = [];

%Threshold for separating gap from break in seconds
boutBreak = 0.25; 

%Extract fields from the annotation elements
segTypes = getAnnotationCell(elements,'segType');
sTime = getAnnotationCell(elements,'segFileStartTimes');
eTime = getAnnotationCell(elements,'segFileEndTimes');

%Unfold cell arrays to vector arrays
for i = 1:length(elements)
    sylLabels = [sylLabels; segTypes{i}];
    starts = [starts; sTime{i}'];
    ends = [ends; eTime{i}'];
    
    %While in the loop, grab gap durations and transition types
    numSylls = length(segTypes{i});
    durs = sTime{i}(2:end)-eTime{i}(1:end-1);
    gapDurs = [gapDurs; durs(durs<=boutBreak)'];
    transLabels = [(segTypes{i}(1:end-1)),(segTypes{i}(2:end))];
    gapLabels = [gapLabels; transLabels(durs<=boutBreak,:)];
end

%All syllable durations
sylDurs = ends-starts;

%Parse syllable durations into types
sylTypes = sort(unique(sylLabels));
for i = 1:length(sylTypes)
   sylDursbyType{i} = sylDurs(sylLabels==sylTypes(i));
end

%Parse gap durations into types
startTypes = sort(unique(gapLabels(:,1)));
endTypes = sort(unique(gapLabels(:,2)));
for i = 1:length(startTypes)
    for j = 1:length(endTypes)
        gapDursbyType{i,j}.durs = gapDurs(gapLabels(:,1)==startTypes(i) & gapLabels(:,2)==endTypes(j));
        gapDursbyType{i,j}.gapType = [num2str(startTypes(i)), '-', num2str(endTypes(j))];
    end
end

function [rhythmMean,rhythmStd] = CalcRhythm(audio)
%Calculate the rhythm of the song based on the pressure amplitude
NFFT = 512;
for i = 1:length(audio)
    tmp = mean(running_windows(audio{i}'.^2,220,44)); %Calculate and smooth amplitude
    xformedAud(i,:) = abs(fft(tmp,NFFT)); %Get spectrum
end

%Take basic stats over the whole set
rhythmMean = mean(xformedAud(:,1:100),1);
rhythmStd = std(xformedAud(:,1:100),1);

function [template_Intermed,target,templatesyllBreaks,templatemotifBreak] = createComplexTemplate(data,fs,buffer,starts)
%Select the average-length rendition to use as the template
target = 1;
renditions = size(data,2);
specPowTS = [];
specLength = [];
specPow = [];

for i=1:renditions;
    %Calculate the total power in the rendition spectrogram
    specPow(i) = sum(sum(data{i}));
    
    %Find the motif edges and calculate motif length
    [~, motifBreak] = getEdges(data{i},starts(:,i),buffer,'robust');
    motifLength(i) = motifBreak(2)-motifBreak(1);
end

%Find the most average length rendition for use as the initial template
aveLength = mean(motifLength);

%Takes the most-average length rendition and uses it alone as initial template
[~, target] = min(abs(motifLength-aveLength));
target = target(1);
template_Initial = data{target(1)};

%Warp all spectrograms to the initial template using Cengiz's method
for i=1:renditions
     %[DTWdist p] = DTWsubWeightedNormBandFast(template_Initial,data{i});
     [DTWdist p] = DTWWeightedNormBandFastest(template_Initial,data{i});
     [audioCube(i,:,:)] = alignSeries(data{i},p);
end

%Take the mean across all time-freq bins
Spec = squeeze(mean(audioCube,1));

%Scale Power to match sample average
template_Intermed = Spec.*(mean(specPow)/sum(Spec(:)));

%Find the syllable and motif edges for the derived template
[templatesyllBreaks,templatemotifBreak] = getEdges(template_Intermed,starts(:,target),buffer,[]);

%Output the data to figure for spot check
figure
scatterhist(1:length(motifLength),motifLength,'NBins',25)
hold on
scatter(target,motifLength(target))
line([1,length(motifLength)],[aveLength,aveLength],'Color','r','LineStyle','--')
title('Initial Motif Length Measurements')
xlabel('Rendition')
ylabel('Motif Length (ms)')
axis tight

% figure
% imagesc(-1*template_Intermed)
% hold on
% for i = 1:length(templatesyllBreaks(:))
%     line([templatesyllBreaks(i), templatesyllBreaks(i)],[0,size(template_Intermed,1)],'Color','k','LineWidth',2)
% end
% axis tight; axis xy;
% title('Derived Template for Alignment')
% xlabel('Time (ms)')
% ylabel('Freq Bins')

function [syllBreaks,motifBreak] = getEdges(data,starts,buffer,type)

Crossings.Up = [];
Crossings.Down = [];
    
%Sum across frequency bins for each rendition
specPowTS = -1*sum(data);
specPow = sum(sum(data));

%Rerun EM until a solution is found
iterations = 1;
gMix = [];

while (isempty(Crossings.Up) || isempty(Crossings.Down))
    gMix = gmdistribution.fit(specPowTS',2);
    [powerM,pnt] = sort(gMix.mu);
    temp = sqrt(squeeze(gMix.Sigma));
    powerStd = temp(pnt);
    thresh = powerM(1)+1*powerStd(1); %1SDs above the silent Gaussian mean

    %Find the crossing points to estimate the onset and offsets of
    %motifs
    ind = 2:length(specPowTS);
    Crossings.Up = find(specPowTS(ind)>thresh & specPowTS(ind-1)<thresh);
    Crossings.Down = find(specPowTS(ind)<thresh & specPowTS(ind-1)>thresh);
    
    if iterations > 5
        print(['gmdistribution.fit has been run ' num2str(iterations) ' without converging on a solution.'])
        return
    end
    iterations = iterations + 1;
end

%Use the starting positions to capture the start/stop points for all
%syllables of the rendition
syllBreaks = [];
motifBreak = [];

onStarts = (starts(1:2:end)-starts(1)+1000*buffer);
offStarts = (starts(2:2:end)-starts(1)+1000*buffer);

 for k = 1:length(onStarts)
     [offset, pntr] = min(abs(Crossings.Up-onStarts(k)));
      if offset > 10 && strcmp(type,'robust')
         syllBreaks(k,1) = onStarts(k);
      else
        syllBreaks(k,1) = Crossings.Up(pntr);
      end
     
     [offset, pntr] = min(abs(Crossings.Down-offStarts(k)));
      if offset > 10 && strcmp(type,'robust')
         syllBreaks(k,2) = offStarts(k);
      else
        syllBreaks(k,2) = Crossings.Down(pntr);
      end
 end
motifBreak(1) = syllBreaks(1,1); motifBreak(2) = syllBreaks(end,2);

function [filtAudio] = Prep(handles,audio)
if ~isempty(audio)
    %Constants for Bandpass Audio (300-10000kHz)
    HP_fNorm = 300/(44150/2);
    LP_fNorm = 6500/(44150/2); %Changed 7/26/12 to match Farhan and Cengiz
    [BP_b,BP_a] = butter(2,[HP_fNorm LP_fNorm]);
    renditions = size(audio,2);
end

parfor i = 1:renditions
    if ~isempty(audio)
        filtAudio{i} = filtfilt(BP_b,BP_a,audio{i});
    end
end

function SylFeatureSet = parseSpectral(SylFeatureSet,rendFeatures,elements,sylTypes)

%Get the list of all fields included in the features structurs
names = fieldnames(rendFeatures);

%Loop through each of the annotated syllables, parse the feature vectors,
%and file in to the appropriate structure
for i = 1:length(elements.segType)
    %Get indices and pointers to pull and push the right data
    ind = floor(1000*elements.segFileStartTimes(i)):ceil(1000*elements.segFileEndTimes(i));
    pntr = find(sylTypes == elements.segType(i),1);
    if isempty(pntr)
        error('Syllable types do not match. Check in parseSpectral function.')
        return
    end
    
    try
        listpntr = length(SylFeatureSet(pntr).ts)+1;
    catch
        listpntr = 1;
    end
    
    %Work through each field to parse the data
    SylFeatureSet(pntr).mean(listpntr).dur= length(ind);
    for j = 1:length(names)
        %For ease of reading, copy rendition-length feature to tmp var
        tmp = eval(['rendFeatures.' names{j}]);
        eval(['SylFeatureSet(pntr).ts{listpntr}.' names{j} ' = tmp(ind);']);
        eval(['SylFeatureSet(pntr).mean(listpntr).' names{j} ' = mean(tmp(ind));']);
    end
    %SylFeatureSet(pntr).dur(listpntr)= length(ind);
    SylFeatureSet(pntr).ID = elements.segType(i);
end

function [X,Ys] = epdf(data,nbins,minXin,maxXin)
%Take in data vector and create probability distribution function given the
%bins and min/max values

xin = reshape(data, numel(data), 1 );
if ~isreal( xin )
    xin = abs( xin );
end

if floor( nbins ) ~= nbins
    error( 'Number of bins should be integer value' );
end
if nbins < 2
    error( 'Number of bins should be positive integer greater than 1 ' );
end

step = (maxXin - minXin) / (nbins-1);
binc = minXin : step : maxXin;     
[N, X] = hist(xin, binc);
Ys = N/sum(N);

function FeatArray = getFeatArray(FeatStruct,sylPntr,field)

numSyls = length(FeatStruct(sylPntr).mean);
FeatArray = [];
for i = 1:numSyls
    FeatArray(i) = getfield(FeatStruct(sylPntr).mean(i),char(field));
end

function updateTimeDistr(handles,hourVect,minVect)

binWidth = 0.5; %width of bins in hours
ToD = hourVect+minVect/60;
hs = histc(ToD,0:binWidth:24);

axes(handles.axes_timeDistr)
h = bar(.5:binWidth:23.5,hs(2:end-1));
xlim([0 24]);
set(h,'facecolor',[.5 .5 .5])
box off;
xlabel('Time of day');
ylabel('File Count')

function sylTypes = updateTypeDistr(handles,segTypes)
%Collapse structure into single array
totTypes = [];
for i = 1:length(segTypes)
    totTypes = [totTypes;segTypes{i}];
end

%Parse array for type counts
sylTypes = sort(unique(totTypes));
hs = histc(totTypes,sylTypes);
for i = 1:length(sylTypes)
    xticklabel{i} = num2str(sylTypes(i));
end

axes(handles.axes_sylDistr);
h = bar(1:length(sylTypes),hs);
set(handles.axes_sylDistr,'XTickLabel',xticklabel,'XTick',1:length(sylTypes))
xlim([0 length(sylTypes)+1]);
set(h,'facecolor',[.5 .5 .5]);
box off;
xlabel('Syllable Type');
ylabel('Count')
hold off

function handles = updatefeatStats(handles)

featNames =  fieldnames(handles.data.spectral.syl_features(1).mean(1));
feat = featNames(get(handles.popup_featStats,'Value'));

%Generate feature arrays for each feature
axes(handles.axes_featStats)
set(handles.axes_featStats,'NextPlot','ReplaceChildren')
cla
hold on
for i = 1:length(handles.data.spectral.syl_featuresKEY)
    output = getFeatArray(handles.data.spectral.syl_features,i,feat);
    stds(i) = std(output);
    means(i) = mean(output);
    maxes(i) = max(output);
    mins(i) = min(output);
    tempTypes{i} = num2str(handles.data.spectral.syl_featuresKEY(i));
    
    jitter = randn(1,length(output))*0.1;
    xs = (ones(1,length(output)).*i)+jitter;
    scatter(xs,output,'.')
end
errorbar(1:length(handles.data.spectral.syl_featuresKEY),means,stds,'ok','LineWidth',1.5)
ytick = linspace(.9*min(mins),1.1*max(maxes),3);
for i=1:length(ytick)
    ylabel{i} = num2str(ytick(i),3);
end
set(handles.axes_featStats,'FontSize',10)
xlim([0,length(handles.data.spectral.syl_featuresKEY)+1])
set(handles.axes_featStats,'XTickLabel',tempTypes,'XTick',1:length(handles.data.spectral.syl_featuresKEY))
set(handles.axes_featStats,'YTickLabel',ylabel,'YTick',ytick)
hold off

function handles = updateFeatVFeat(handles)

featNames =  fieldnames(handles.data.spectral.syl_features(1).mean(1));
feat1 = featNames(get(handles.popup_feat1,'Value'));
feat2 = featNames(get(handles.popup_feat2,'Value'));

%Generate feature arrays for each feature
axes(handles.axes_featVfeat)
set(handles.axes_featVfeat,'NextPlot','ReplaceChildren')
cla
hold on
for i = 1:length(handles.data.spectral.syl_featuresKEY)
    output1 = getFeatArray(handles.data.spectral.syl_features,i,feat1);
    maxes1(i) = max(output1);
    mins1(i) = min(output1);
    
    output2 = getFeatArray(handles.data.spectral.syl_features,i,feat2);
    maxes2(i) = max(output2);
    mins2(i) = min(output2);
    
    scatter(output1,output2,'.')
end

xlim([min(mins1)-.1*min(mins1),1.1*max(maxes1)])
ylim([min(mins2)-.1*min(mins2),1.1*max(maxes2)])

xtick = linspace(min(mins1),max(maxes1),3);
for i=1:length(xtick)
    xlabel{i} = num2str(xtick(i),3);
end

ytick = linspace(min(mins2),max(maxes2),3);
for i=1:length(ytick)
    ylabel{i} = num2str(ytick(i),3);
end

set(handles.axes_featVfeat,'FontSize',10)
set(handles.axes_featVfeat,'XTick',xtick,'XTickLabel',xlabel)
set(handles.axes_featVfeat,'YTick',ytick,'YTickLabel',ylabel)

hold off

function handles = updateSimStats(handles)
%This is the number of syllables to be randomly chosen for the "quick"
%analysis
quickNum = 20;

%The amount of jitter/slack (in seconds) allowable in MatchScore
slack = 0.01;

%Generate feature arrays for each feature
axes(handles.axes_simStats)
set(handles.axes_simStats,'NextPlot','ReplaceChildren')
cla
hold on
handles.data.spectral.matchScore = [];
for i = 1:length(handles.data.spectral.syl_featuresKEY)
    if get(handles.check_quickSim,'Value') && length(handles.data.spectral.syl_features(i).ts)>quickNum
        index = sort(randperm(length(handles.data.spectral.syl_features(i).ts),quickNum));
    else
        index = 1:length(handles.data.spectral.syl_features(i).ts);
    end
    ele = 0;
    for j = 1:length(index)
        for k = 1:j
            %Save as an matrix to preserve pairwise comparisons
            handles.data.spectral.matchScore{i}.matrix(j,k) = sylMatchScore(handles.data.spectral.syl_features(i).ts{index(j)}, handles.data.spectral.syl_features(i).ts{index(k)}, slack, false);
            %Save as an array (without useless values) for ease of coding
            if j~=k
                ele= ele+1;
                handles.data.spectral.matchScore{i}.array(ele) = handles.data.spectral.matchScore{i}.matrix(j,k);
            end
        end
    end
    means(i) = mean(handles.data.spectral.matchScore{i}.array);
    stds(i) = std(handles.data.spectral.matchScore{i}.array);
    
    tempTypes{i} = num2str(handles.data.spectral.syl_featuresKEY(i));
    
    jitter = randn(1,ele)*0.1;
    xs = (ones(1,ele).*i)+jitter;
    scatter(xs,handles.data.spectral.matchScore{i}.array,'.')
end
sylInds = find(handles.data.spectral.syl_featuresKEY<100);
%errorbar(1:(length(handles.data.spectral.syl_featuresKEY)+1),[means,mean(means(sylInds))],[stds,std(means(sylInds))],'ok','LineWidth',1.5)
errorbar(1:length(handles.data.spectral.syl_featuresKEY),means,stds,'ok','LineWidth',1.5)
errorbar(length(handles.data.spectral.syl_featuresKEY)+1,mean(means(sylInds)),std(means(sylInds)),'sr','LineWidth',1.5)
ytick = linspace(0,1,5);
for i=1:length(ytick)
    ylabel{i} = num2str(ytick(i),2);
end
tempTypes{end+1} = 'M';

set(handles.axes_simStats,'FontSize',10)
xlim([0,length(handles.data.spectral.syl_featuresKEY)+1])
ylim([-.1,1.1])
set(handles.axes_simStats,'XTickLabel',tempTypes,'XTick',1:(length(handles.data.spectral.syl_featuresKEY)+1))
set(handles.axes_simStats,'YTickLabel',ylabel,'YTick',ytick)

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_addFolder_Callback(hObject, eventdata, handles)
% hObject    handle to push_addFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get user input for where to find the .dat files and check location
directory_name = uigetdir('Choose location of recordings to analyze:');

if isequal(directory_name,0) || ~isstr(directory_name)
    set(handles.text_message,'String','Data location invalid or cancelled. Pick a valid folder.');
    guidata(hObject, handles);
    return
end

cd(directory_name)
filelistDAT = dir('*.dat');
filelistWAV = dir('*.wav');
filelist = [filelistDAT;filelistWAV];

if isempty(filelist)
    set(handles.text_message,'String',[directory_name ' does not contain any *.dat or *.wav files. Pick a valid folder.']);
    guidata(hObject, handles);
    return
else
    %If they don't hit cancel, update data location...

    %This is the code for handling multiple cell files. 
    if ~isempty(handles.dirname)
        button = questdlg('Do you want to add this data location to the previously loaded one(s) or replace them?','Add a data location','Add to current','Replace current','Add to current');
        if (strcmp(button,'Replace current'))
            %Overwrite the current directory handle
            handles.dirname = directory_name;
            
            %Add directory location to the file structure
            for i = 1:length(filelist)
                filelist(i).dirname = directory_name;
            end
            handles.totFilelist = filelist;
            
            %Reset all the tags and objects that should go back to initial values 
            %Set the message to show dataset expanded
            set(handles.text_message,'String',['Files from from ' directory_name ' are the new dataset.']);    
        elseif (strcmp(button,'Add to current'))
            if strcmp(directory_name,handles.dirname)
                set(handles.text_message,'String',[directory_name ' is already included in the dataset.']);
                guidata(hObject, handles);
                return
            else
                %handles.dirname = {handles.dirname; directory_name};
                
                if size(handles.dirname,1)==1
                    handles.dirname = {handles.dirname;directory_name};
                else
                    handles.dirname{length(handles.dirname)+1} = directory_name;
                end
            end
            
            %Add directory location to the file structure
            for i = 1:length(filelist)
                filelist(i).dirname = directory_name;
            end
            handles.totFilelist =  [handles.totFilelist; filelist];
            
            %Set the message to show dataset expanded
            set(handles.text_message,'String',['Files from from ' directory_name ' added to the dataset.']);
        end
    %This is the code for loading the first/only folder  
    else
        %Overwrite the current directory handle
        handles.dirname = directory_name;
        
        %Add directory location to the file structure
        for i = 1:length(filelist)
            filelist(i).dirname = directory_name;
        end
        handles.totFilelist = filelist;

        %Reset all the tags and objects that should go back to initial values 

        %Set the message to show dataset expanded
        set(handles.text_message,'String',['Files from from ' directory_name ' are the new dataset.']);
    end
    set(handles.listbox_data,'String',handles.dirname);
end

%If the annotation file is already loaded, cross reference the data
if ~isempty(handles.elements)
    %Pare down the totFilelist set to only those files that have a record
    %in the loaded annotations; dump the results in handles.filelist
    [handles.filelist] = XrefAnnotData(handles.totFilelist, handles.keys);
    set(handles.listbox_keys,'String', handles.keys)
    for i=1:length(handles.filelist)
        names{i} = handles.filelist(i).name;
    end
    set(handles.listbox_filelist,'String',names)
        
    %Update the record selection criteria
    %.........the options for sorting by drug status
    contents = ['All' sort(unique(handles.drugstatus))];
    set(handles.popup_drugstatus,'String',contents);
    set(handles.popup_drugstatus,'Value',1);
    handles.filter_DrugStat = contents(1);

    %.........the options for sorting by directed status
    contents = ['All' sort(unique(handles.directstatus))];
    set(handles.popup_directstatus,'String',contents);
    set(handles.popup_directstatus,'Value',1);
    handles.filter_DirectStat = contents(1);

    %.........the upper and lower bounds of the filenum selector
    set(handles.edit_fromRec,'String',num2str(min(handles.filenums)));
    handles.filter_fromRec = min(handles.filenums);
    set(handles.edit_toRec,'String',num2str(max(handles.filenums)));
    handles.filter_toRec = max(handles.filenums);

    %.........the upper and lower bounds of the time selector
    sTime = [num2str(handles.hourStamp(1)), num2str(handles.minuteStamp(1))];
    eTime = [num2str(handles.hourStamp(end)), num2str(handles.minuteStamp(end))];
    set(handles.edit_fromTime,'String',sTime);
    handles.filter_fromTime = [handles.hourStamp(1), handles.minuteStamp(1)];
    set(handles.edit_toTime,'String',eTime);
    handles.filter_toTime = [handles.hourStamp(end), handles.minuteStamp(end)];
end

guidata(hObject, handles);

function push_addAnnot_Callback(hObject, eventdata, handles)

[file,path] = uigetfile('*.mat', 'Choose an audio annotation .mat file to load:');

if isequal(file,0) || isequal(path,0)
    set(handles.text_message,'String','Annotation location invalid or cancelled. Pick a valid file.');
else
    %if they don't hit cancel, update annotation location and load the file
    if ~isempty(handles.annotationFilename)
        button = questdlg('Do you want to add this annotation to the current list or replace them?','Add an annotation','Add to current','Replace current','Add to current');
        if (strcmp(button,'Replace current'))
            handles.annotationFilename = [path,file];
            handles.annotitles = file;
            annotation = aaLoadHashtable(handles.annotationFilename);

            %Copy out the two main structures of the annotation file
            handles.elements = annotation.elements;
            handles.keys = annotation.keys;
            
            %Sort the keys alpha-numerically (which is also temporally) and
            %use the index to sort the elements vector as well
            [handles.keys ind] = sort(handles.keys);
            handles.elements = handles.elements(ind);
            
            %Parse out the sortable fields of the elements structure
            handles.filenums=getAnnotationVector(handles.elements,'filenum');
            handles.drugstatus=getAnnotationVector(handles.elements,'drugstatus');
            handles.directstatus=getAnnotationVector(handles.elements,'directstatus');
            [~,handles.hourStamp] = getKeysSubset(handles.keys,'h');
            [~,handles.minuteStamp] = getKeysSubset(handles.keys,'min');
            
            %Set the message to show annotation replaced
            set(handles.listbox_annotation,'String',handles.annotitles)
            set(handles.text_message,'String',[file ' is the new annotation.']);
            
        elseif (strcmp(button,'Add to current'))
            if strcmp(file,handles.annotitles)
                set(handles.text_message,'String',[file ' is already included in the dataset.']);
                guidata(hObject, handles);
                return
            end
            
            if size(handles.annotitles,1)==1
                handles.annotationFilename = {handles.annotationFilename; [path,file]};
                handles.annotitles = {handles.annotitles;file};
            else
                handles.annotationFilename{length(handles.annotationFilename)+1} = [path,file];
                handles.annotitles{length(handles.annotitles)+1} = file;
            end
            
            annotation = aaLoadHashtable([path,file]);
            
            %Copy out the two main structures of the annotation file
            handles.elements = [handles.elements annotation.elements];
            handles.keys = [handles.keys annotation.keys];
            
            %Sort the keys alpha-numerically (which is also temporally) and
            %use the index to sort the elements vector as well
            [handles.keys ind] = sort(handles.keys);
            handles.elements = handles.elements(ind);

            %Parse out the sortable fields of the elements structure
            handles.filenums = getAnnotationVector(handles.elements,'filenum');
            handles.drugstatus = getAnnotationVector(handles.elements,'drugstatus');
            handles.directstatus = getAnnotationVector(handles.elements,'directstatus');
            [~,handles.hourStamp] = getKeysSubset(handles.keys,'h');
            [~,handles.minuteStamp] = getKeysSubset(handles.keys,'min');
            
            %Set the message to show annotation replaced
            set(handles.listbox_annotation,'String',handles.annotitles)
            set(handles.text_message,'String',[file ' is added to the annotation list.']);
        end
    %This is the code for loading the first/only folder    
    else
        handles.annotationFilename = [path,file];
        handles.annotitles = file;
        annotation = aaLoadHashtable(handles.annotationFilename);

        %Copy out the two main structures of the annotation file
        handles.elements = annotation.elements;
        handles.keys = annotation.keys;

        %Sort the keys alphabetically (which is also temporally) and
        %use the index to sort the elements vector as well
        [handles.keys ind] = sort(handles.keys);
        handles.elements = handles.elements(ind);
            
        %Parse out the sortable fields of the elements structure
        handles.filenums=getAnnotationVector(handles.elements,'filenum');
        handles.drugstatus=getAnnotationVector(handles.elements,'drugstatus');
        handles.directstatus=getAnnotationVector(handles.elements,'directstatus');
        [~,handles.hourStamp] = getKeysSubset(handles.keys,'h');
        [~,handles.minuteStamp] = getKeysSubset(handles.keys,'min');
            
        %Set the message to show annotation replaced
        set(handles.listbox_annotation,'String',handles.annotitles)
        set(handles.text_message,'String',[handles.annotationFilename ' is the new annotation.']);
    end
    
    %If a data location is already set, cross-reference the two
    if ~isempty(handles.totFilelist)
        %Pare down the totFilelist set to only those files that have a record
        %in the loaded annotations; dump the results in handles.filelist
        [handles.filelist] = XrefAnnotData(handles.totFilelist, handles.keys);
        set(handles.listbox_keys,'String', handles.keys)
        for i=1:length(handles.filelist)
            names{i} = handles.filelist(i).name;
        end
        set(handles.listbox_filelist,'String',names)

        %Update the record selection criteria
        %.........the options for sorting by drug status
        contents = ['All' sort(unique(handles.drugstatus))];
        set(handles.popup_drugstatus,'String',contents);
        set(handles.popup_drugstatus,'Value',1);
        handles.filter_DrugStat = contents(1);
        
        %.........the options for sorting by directed status
        contents = ['All' sort(unique(handles.directstatus))];
        set(handles.popup_directstatus,'String',contents);
        set(handles.popup_directstatus,'Value',1);
        handles.filter_DirectStat = contents(1);
        
        %.........the upper and lower bounds of the filenum selector
        set(handles.edit_fromRec,'String',num2str(min(handles.filenums)));
        handles.filter_fromRec = min(handles.filenums);
        set(handles.edit_toRec,'String',num2str(max(handles.filenums)));
        handles.filter_toRec = max(handles.filenums);
        
        %.........the upper and lower bounds of the time selector
        sTime = [num2str(handles.hourStamp(1),'%02d'), num2str(handles.minuteStamp(1),'%02d')];
        eTime = [num2str(handles.hourStamp(end),'%02d'), num2str(handles.minuteStamp(end),'%02d')];
        set(handles.edit_fromTime,'String',sTime);
        handles.filter_fromTime = [handles.hourStamp(1), handles.minuteStamp(1)];
        set(handles.edit_toTime,'String',eTime);
        handles.filter_toTime = [handles.hourStamp(end), handles.minuteStamp(end)];
    end
end
guidata(hObject, handles);

function push_clearAnnots_Callback(hObject, eventdata, handles)

button = questdlg('Are you sure you want to clear all annotations?','Clear annotations?','Clear''em All!','Nooooo!','Nooooo!');

if (strcmp(button,'Clear''em All!'))
    handles.annotationFilename = [];
    handles.annotitles = [];
    handles.filelist = [];
    handles.dataset = [];

    %Clear the two main structures of the annotation file
    handles.elements = [];
    handles.keys = [];

    %Parse out the sortable fields of the elements structure
    handles.filenums=[];
    handles.drugstatus=[];
    handles.directstatus=[];

    %Set the message to show annotation replaced
    set(handles.listbox_annotation,'String',handles.annotitles)
    set(handles.text_message,'String','All annotations cleared.');
else
    set(handles.text_message,'String','Annotations are unchanged.');
end

guidata(hObject, handles);

function push_clearFolders_Callback(hObject, eventdata, handles)
button = questdlg('Are you sure you want to clear all folders?','Clear folder?','Clear''em All!','Nooooo!','Nooooo!');

if (strcmp(button,'Clear''em All!'))
    handles.dirname = [];

    %Clear the two main repositories of the file info
    handles.totFilelist = [];
    handles.filelist = [];
    handles.dataset = [];
    
    %Set the message to show annotation replaced
    set(handles.listbox_data,'String',handles.dirname);
    set(handles.text_message,'String','All folders cleared.');
else
    set(handles.text_message,'String','Data folders are unchanged.');
end

guidata(hObject, handles);

function push_loadSylIndex_Callback(hObject, eventdata, handles)
[file,dpath]=uigetfile('*.mat');

if isequal(file,0) || isequal(dpath,0)
    set(handles.text_message,'String','Operation canceled or selected file is invalid.');
else
    load([dpath file]);
    handles.templates=templates;
    set(handles.text_message,'String',['Syllable Index was loaded from ' file '.']);
end
guidata(hObject, handles);

function push_createDataset_Callback(hObject, eventdata, handles)
%Load the data from the saved .dat and .wav files to the workspace
if isempty(handles.filelist)
    set(handles.text_message,'String','No files found to load. Check annotations and data folders.')
else
    %Work through each file to extract the audio data
    handles.dataset = [];
    for i=1:length(handles.filelist)
        if ~strcmp(handles.filelist(i).name,'DUMMY') %If the record isn't crap...
            if strcmp(handles.filelist(i).name(end-2:end),'dat') %...and it's a multi-channel dat file
                [channels, fs] = getChannels([handles.filelist(i).dirname filesep handles.filelist(i).name]);
                handles.dataset{i} = channels(1,:); %Take the audio channel only
            elseif strcmp(handles.filelist(i).name(end-2:end),'wav') %...and it's a single channel wav file.
                [handles.dataset{i}, fs] = wavread([handles.filelist(i).dirname filesep handles.filelist(i).name]);
            end
        end
    end
    set(handles.text_message,'String','Dataset creation now completed')
end

if get(handles.check_useAllSyls,'Value')
    %Clear out all of the plots and data
    handles = clearAllDataAxes(handles);
    handles.data = [];

    %Create a record index based on user filter selections
    [handles.filtInd] = filterRecords(handles);

    %Update the Annotation Stats Plots
    updateTimeDistr(handles,handles.hourStamp(handles.filtInd==1),handles.minuteStamp(handles.filtInd==1))
    segTypes=getAnnotationCell(handles.elements,'segType');
    handles.data.sylTypes = updateTypeDistr(handles,segTypes(handles.filtInd==1));
    
    for i = 1:length(handles.data.sylTypes)
        contents{i} = num2str(handles.data.sylTypes(i));
    end
    set(handles.popup_sylSelect,'String',contents)
end

guidata(hObject, handles);

function listbox_annotation_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox_data_Callback(hObject, eventdata, handles)

function listbox_data_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox_keys_Callback(hObject, eventdata, handles)

function listbox_keys_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox_filelist_Callback(hObject, eventdata, handles)

function listbox_filelist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Record Selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_fromRec_Callback(hObject, eventdata, handles)

if str2double(get(handles.edit_fromRec,'String')) >= min(handles.filenums) && str2double(get(handles.edit_fromRec,'String')) <= str2double(get(handles.edit_toRec,'String'))
    handles.filter_fromRec = str2double(get(handles.edit_fromRec,'String'));
else
    set(handles.edit_fromRec,'String',num2str(min(handles.filenums)));
    handles.filter_fromRec = min(handles.filenums);
end
guidata(hObject, handles);

function edit_fromRec_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_toRec_Callback(hObject, eventdata, handles)

if str2double(get(handles.edit_toRec,'String')) <= max(handles.filenums) && str2double(get(handles.edit_toRec,'String')) >= str2double(get(handles.edit_fromRec,'String'))
    handles.filter_toRec = str2double(get(handles.edit_toRec,'String'));
else
    set(handles.edit_toRec,'String',num2str(max(handles.filenums)));
    handles.filter_toRec = max(handles.filenums);
end
guidata(hObject, handles);

function edit_toRec_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_directstatus_Callback(hObject, eventdata, handles)

contents = get(handles.popup_directstatus,'String');
handles.filter_DirectStat = contents(get(handles.popup_directstatus,'Value'));

guidata(hObject, handles);

function popup_directstatus_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_drugstatus_Callback(hObject, eventdata, handles)

contents = get(handles.popup_drugstatus,'String');
handles.filter_DrugStat = contents(get(handles.popup_drugstatus,'Value'));

guidata(hObject, handles);

function popup_drugstatus_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_fromTime_Callback(hObject, eventdata, handles)

function edit_fromTime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_toTime_Callback(hObject, eventdata, handles)

function edit_toTime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Motif Selections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function popup_seqSylls_Callback(hObject, eventdata, handles)

%Returns the numerical sequence that the user has selected
chosenSeq = handles.seq(get(handles.popup_seqSylls, 'Value')); %chosen sequence as number
handles.Seq = chosenSeq;

%Finds the index of the chosen sequence as they appear in the
%handles.all_seq catalog
startIndex = find(handles.allSeq==chosenSeq);

%Gives the filtered record/annotation number and the syllable number within
%the annotation
handles.chosenStartSeq = [];
for i=1:length(startIndex)
    handles.chosenStartSeq(i,:) = handles.startSeq(startIndex(i),:); 
end

%Filter the elements array 
filt_elements = handles.elements(handles.filtInd == 1);

%Determine the selected sequence and prep it for displaying the appropriate
%template on axes_template
segType = filt_elements{handles.chosenStartSeq(1,1)}.segType;
handles.sequence=[];
handles.sequence(1:handles.numSylls) = segType(handles.chosenStartSeq(1,2):handles.chosenStartSeq(1,2)+handles.numSylls-1);

%Clear out old data and templates related to the previous selections
handles.data = [];

guidata(hObject, handles);

function popup_seqSylls_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_numSylls_Callback(hObject, eventdata, handles)

if (~isfield(handles,'filter_fromRec') || ~isfield(handles,'filter_toRec') || ~isfield(handles,'filter_DrugStat') || ~isfield(handles,'filter_DirectStat'))
    warndlg('Select the records that you want to display before continuing.');
    uiwait;
elseif isempty(handles.dataset)
    set(handles.text_message,'String','There is no dataset in memory.  Create dataset before continuing.')
    set(handles.popup_numSylls,'Value',1);
else
    %Copy out the user's selection from the popup menu
    handles.numSylls = get(handles.popup_numSylls, 'Value');

    %Clear out all of the plotted data
    handles = clearAllDataAxes(handles);
    handles.data = [];
    
    %Uncheck 'Use All' box
    set(handles.check_useAllSyls,'Value',0);
    
    %Create a record index based on user filter selections
    [handles.filtInd] = filterRecords(handles);
    
    %Run through the now refined set of records to extract the sequences
    [handles.allSeq handles.startSeq]=sortSequences(handles, handles.numSylls); %sort the sequences of handles.numSyll in terms of how often they appear

    %Sort and display sequences by prevalance
    histx=(1:10^handles.numSylls);
    hist_seq=hist(handles.allSeq,histx);
    
    %Choose the 10 most abundant sequences to display in the sequence pop-up window
    [sortSeq, sortIndex]=sort(hist_seq,'descend');
    numSeq=length(find(sortSeq>0));
    if isempty(numSeq)
        warndlg('No such sequences found.');
        uiwait;
        return;
    elseif (numSeq<10)
        popupNum=numSeq;
    else
        popupNum=10;
    end
    popupSeqString(1,:)=[num2str(sortIndex(1)),' (n=', ' ',num2str(sortSeq(1)),')'];
    
    handles.seq=sortIndex(1);
    if (popupNum>1)
        for i=2:popupNum %display the frequency of the various sequences)
            padLength=length(num2str(sortSeq(1)))-length(num2str(sortSeq(i)));
            if (padLength==1)
                freq(i,:)=['  ',num2str(sortSeq(i))];
            elseif (padLength==2)
                freq(i,:)=['   ',num2str(sortSeq(i))];
            else
                freq(i,:)=[' ',num2str(sortSeq(i))];
            end
            handles.seq(i)=sortIndex(i);
            popupSeqString(i,:)=[num2str(sortIndex(i)),' (n=', freq(i,:),')'];
        end
    end
        
    set(handles.popup_seqSylls,'String',popupSeqString);
    set(handles.popup_seqSylls,'Value',1);
end

guidata(hObject, handles);

function popup_numSylls_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function check_useAllSyls_Callback(hObject, eventdata, handles)
if isempty(handles.dataset)
    set(handles.text_message,'String','There is no dataset in memory.  Create dataset before continuing.')
    set(handles.popup_numSylls,'Value',1);
else
    if ~get(handles.check_useAllSyls,'Value') %If user unchecks the box, reset the motif pop-ups
        %Clear out all of the plots and data
        handles = clearAllDataAxes(handles);
        handles.data = [];
        
    else %Otherwise, run through the dataset filtering process...
        %Clear out all of the plots and data
        handles = clearAllDataAxes(handles);
        handles.data = [];
    
        %Create a record index based on user filter selections
        [handles.filtInd] = filterRecords(handles);
        
        %Update the Annotation Stats Plots
        updateTimeDistr(handles,handles.hourStamp(handles.filtInd==1),handles.minuteStamp(handles.filtInd==1))
        segTypes=getAnnotationCell(handles.elements,'segType');
        handles.data.sylTypes = updateTypeDistr(handles,segTypes(handles.filtInd==1));
        
        for i = 1:length(handles.data.sylTypes)
            contents{i} = num2str(handles.data.sylTypes(i));
        end
        set(handles.popup_sylSelect,'String',contents)
    end
end

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectral Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_spectral1_Callback(hObject, eventdata, handles)

%Generate Workset given current selections
[handles.curAudio handles.curElements] = getWorkset(handles);

%Calculate Spectral Features for each rendition; parse by syllable
SylFeatureSet = [];
for i = 1:length(handles.curAudio)
    handles.data.spectral.rend_features{i} = koenigSpectral(handles.curAudio{i},handles.fs);
    SylFeatureSet = parseSpectral(SylFeatureSet,handles.data.spectral.rend_features{i},handles.curElements{i},handles.data.sylTypes);
end
handles.data.spectral.syl_featuresKEY = handles.data.sylTypes;
handles.data.spectral.syl_features = SylFeatureSet;

%Update featStat plot given the user selections
handles = updatefeatStats(handles);

%Update featVfeat plot given the user selections
handles = updateFeatVFeat(handles);

guidata(hObject, handles);

function push_spectral2_Callback(hObject, eventdata, handles)

%Update featVfeat plot given the user selections
handles = updateSimStats(handles);

guidata(hObject, handles);

function popup_featStats_Callback(hObject, eventdata, handles)

%Update featStat plot given the user selections
updatefeatStats(handles)

guidata(hObject, handles);

function popup_featStats_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_feat1_Callback(hObject, eventdata, handles)
%Update featVfeat plot given the user selections
handles = updateFeatVFeat(handles);

guidata(hObject, handles);

function popup_feat1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_feat2_Callback(hObject, eventdata, handles)
%Update featVfeat plot given the user selections
handles = updateFeatVFeat(handles);

guidata(hObject, handles);

function popup_feat2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function check_quickSim_Callback(hObject, eventdata, handles)
%Checking this box, sets the paired box to false
set(handles.check_fullSim,'Value',0);

%You may not uncheck the box
if ~get(handles.check_quickSim,'Value')
    set(handles.check_quickSim,'Value',1);
end
guidata(hObject, handles);

function check_fullSim_Callback(hObject, eventdata, handles)
%Checking this box, sets the paired box to false
set(handles.check_quickSim,'Value',0);

%You may not uncheck the box
if ~get(handles.check_fullSim,'Value')
    set(handles.check_fullSim,'Value',1);
end
guidata(hObject, handles);

function push_createSylAve_Callback(hObject, eventdata, handles)

%Clear the axes
axes(handles.axes13)
hold off
cla

%Determine which Syllable has been selected.
sylTypes = get(handles.popup_sylSelect,'String');
selectedSyl = str2num(char(sylTypes(get(handles.popup_sylSelect,'Value'))));
ind = get(handles.popup_sylSelect,'Value');
handles.data.sylSnips{ind} = [];

%Work through the dataset and extract the audio snippets that correspond to
%this syllable

%Create a record index based on user filter selections
[handles.filtInd] = filterRecords(handles);

%Run through the now refined set of records to extract the sequences
[allSeq startSeq ~]=getSequences(handles,1);

%Finds the index of the chosen sequence as they appear in the catalog
startIndex = find(allSeq==selectedSyl);

%Gives the filtered record/annotation number and the syllable number within
%the annotation
chosenStartSeq = [];
for i=1:length(startIndex)
    chosenStartSeq(i,:) = startSeq(startIndex(i),:); 
end

%Parse data from the dataset into two matrices (audio and neuro)
%Filter the annotation and dataset
filt_elements = handles.elements(handles.filtInd==1);
filt_dataset = handles.dataset(handles.filtInd==1);
buffer = 0.025;
renditions = size(chosenStartSeq,1);
for i=1:renditions
    %Grab the start and end times of the syllable
    handles.data.sylSnips{ind}.start(i) = filt_elements{chosenStartSeq(i,1)}.segFileStartTimes(chosenStartSeq(i,2));
    handles.data.sylSnips{ind}.end(i) = filt_elements{chosenStartSeq(i,1)}.segFileEndTimes(chosenStartSeq(i,2));

    startT = floor(((handles.data.sylSnips{ind}.start(i)-buffer)*handles.fs));
    startT = max(1,startT);
    endT = ceil(((handles.data.sylSnips{ind}.end(i)+buffer)*handles.fs));
    endT = min(endT,length(filt_dataset{chosenStartSeq(i,1)}));
    
    %Grab the segment of data that corresponds to each selected motif
    handles.data.sylSnips{ind}.ts{i} = filt_dataset{chosenStartSeq(i,1)}(startT:endT);
end

%Prep/filter audio and neuro data
PPaudio = Prep(handles,handles.data.sylSnips{ind}.ts);

%Create log Spectrogram matrix for each motif rendition
strt = 4; stp = 94; %These bins correspond to 300-8000Hz
for i=1:renditions
    % Calculate STFT features (5ms window, 1ms advance)
    [S,F,T,P] = spectrogram((PPaudio{i}/(sqrt(mean(PPaudio{i}.^2)))),220,220-44,512,44150);
    handles.data.sylSnips{ind}.rawSpecs{i} = -1*log10(abs(P(strt:stp,:)));
    starts(:,i) = [25,size(P,2)-25];
end

%Generate template for the alignment
% if handles.comTempON
%     %Copy the common template to the active structure
%     handles.data.template = handles.comTemplate;
%     handles.data.templatesyllBreaks = handles.comTemplatesyllBreaks;
%     handles.data.templatemotifBreaks = handles.comTemplatemotifBreaks;
% else
    %Create template for alignment
    [handles.data.sylSnips{ind}.template,~,handles.data.sylSnips{ind}.templatesyllBreaks,~] = createComplexTemplate(handles.data.sylSnips{ind}.rawSpecs,handles.fs,buffer,starts);
%end

%Plot the spectrogram of the template in separate figure
figure
imagesc(-1*handles.data.sylSnips{ind}.template)
hold on
for i = 1:length(handles.data.sylSnips{ind}.templatesyllBreaks(:))
    line([handles.data.sylSnips{ind}.templatesyllBreaks(i), handles.data.sylSnips{ind}.templatesyllBreaks(i)],[0,size(handles.data.sylSnips{ind}.template,1)],'Color','k','LineWidth',2)
end
axis xy; axis tight;
xlabel('Time (ms)')
ylabel('Freq Bins')
title('Derived Template and Syllable Bounds')

%LocLin align all renditions to derived template
[~,templateLength] = size(handles.data.sylSnips{ind}.template);
for i=1:renditions
     %Generate warp path for the audio spectrogram to the template
     [~,rendLength] = size(handles.data.sylSnips{ind}.rawSpecs{i});
     [DTWdist p] = DTWWeightedNormBandFastest(handles.data.sylSnips{ind}.template,handles.data.sylSnips{ind}.rawSpecs{i});
     
     %Get the syllable start/stops based on warping path
     handles.data.sylSnips{ind}.syllBreaks(i,:,:) = getWarpedStarts(p,handles.data.sylSnips{ind}.templatesyllBreaks);
     
     %Syl Onsets & Offsets
     tempT = [];
     tempR = [];
     for j = 1:size(handles.data.sylSnips{ind}.templatesyllBreaks,1)
         tempT = [tempT, handles.data.sylSnips{ind}.templatesyllBreaks(j,:)];
         tempR = [tempR, squeeze(handles.data.sylSnips{ind}.syllBreaks(i,j,:))'];
     end
     tempAnchors = [1, tempT, templateLength];
     rendAnchors = [1, tempR, rendLength];
     
     %Create linear path between the chosen anchor points
     handles.data.sylSnips{ind}.path{i} = [tempAnchors', rendAnchors'];
     
     %Warp the audio spectrograms with the new linear path
     handles.data.sylSnips{ind}.alignedSpecs(i,:,:) = alignSeriesSTW(handles.data.sylSnips{ind}.rawSpecs{i},handles.data.sylSnips{ind}.path{i});
end

%Plot on the GUI axis
axes(handles.axes13)
im = imagesc(-1*squeeze(mean(handles.data.sylSnips{ind}.alignedSpecs,1)));
hold on
for i = 1:length(handles.data.sylSnips{ind}.templatesyllBreaks(:))
    line([handles.data.sylSnips{ind}.templatesyllBreaks(i), handles.data.sylSnips{ind}.templatesyllBreaks(i)],[0,size(handles.data.sylSnips{ind}.template,1)],'Color','k','LineWidth',2)
end
axis xy; axis tight;
xlabel('Time (ms)')
ylabel('Freq Bins')

%Set the flags for pb sensitvity
set(im, 'ButtonDownFcn', @cb_specgram_click);
set(im,'HitTest','on');

%Clear out the old window edges
set(handles.edit_winEdge,'String','')

guidata(hObject, handles);

function push_sylWindow_Callback(hObject, eventdata, handles)
%Find which syllable the popup is set to
ind = get(handles.popup_sylSelect,'Value');

%If present, clear the current window
if isfield(handles.data.sylSnips{ind},'rectHandles')
    if ~isempty(handles.data.sylSnips{ind}.rectHandles)
        delete(handles.data.sylSnips{ind}.rectHandles)
    end
end

if isfield(handles.data.sylSnips{ind},'alignedSpecs') %if the field exists...
    if ~isempty(handles.data.sylSnips{ind}.alignedSpecs) %...and it's not empty for the current syllable
        handles.bWaitingForAddClick = true;
        handles.bWaitingForSeparationClick = false;
        set(handles.push_sylWindow, 'BackgroundColor', 'red');
    else
        set(handles.text_message,'String','Must create an average syllable before setting a window.')
    end
else
    set(handles.text_message,'String','Must create an average syllable before setting a window.')
end

guidata(hObject, handles);

function cb_specgram_click(hObject, evnt)
handles = guidata(hObject);
%fs = handles.fs;
%filenum = handles.filenum;
%mouseMode = get(get(hObject,'Parent'), 'SelectionType');
clickLocation = get(handles.axes13, 'CurrentPoint');

axes(handles.axes13);
if(handles.bWaitingForAddClick)
    %Reset the Window Button
    handles.bWaitingForAddClick = false;
    set(handles.push_sylWindow, 'BackgroundColor', 'white');
    
    %Get the selection region
    rect = rbbox;
    endPoint = get(gca,'CurrentPoint'); 
    point1 = round(clickLocation(1,1:2));       % extract x and y
    point2 = round(endPoint(1,1:2));
    if point1(1)>point2(1)
        point1 = round(endPoint(1,1:2));
        point2 = round(clickLocation(1,1:2));
    end
        
    %Set window edges
    ind = get(handles.popup_sylSelect,'Value');
    handles.data.sylSnips{ind}.winEdge = [point1(1),point2(1)];
    set(handles.edit_winEdge,'String',[num2str(point1(1)) ',' num2str(point2(1))]);
    
    %Draw the window
    [~,c] = size(handles.data.sylSnips{ind}.template);
    handles.data.sylSnips{ind}.rectHandles = patch([point1(1),point2(1),point2(1),point1(1)],[0,0,c,c], [1,1,1], 'EdgeColor', [1,1,1], 'FaceAlpha', 0);
end

guidata(hObject, handles);

function edit_winEdge_Callback(hObject, eventdata, handles)
%Get the input values
temp = str2double(regexp(get(handles.edit_winEdge,'String'),',','split'));

%If the input isn't formatted correctly, kick it back
if length(temp)~=2
    set(handles.text_message,'String','Window must be formatted as: t1,t2');
    return
end
if temp(1)>temp(2)
    set(handles.text_message,'String','Window bound t1 must be less than t2');
    return
end
%Find which syllable the popup is set to
ind = get(handles.popup_sylSelect,'Value');

%If present, clear the current window
if isfield(handles.data.sylSnips{ind},'rectHandles')
    if ~isempty(handles.data.sylSnips{ind}.rectHandles)
        delete(handles.data.sylSnips{ind}.rectHandles)
        handles.data.sylSnips{ind}.rectHandles = [];
    end
end

%Set window edges and draw on the axes
handles.data.sylSnips{ind}.winEdge = [temp(1),temp(2)];
    
[~,c] = size(handles.data.sylSnips{ind}.template);
handles.data.sylSnips{ind}.rectHandles = patch([temp(1),temp(2),temp(2),temp(1)],[0,0,c,c], [1,1,1], 'EdgeColor', [1,1,1], 'FaceAlpha', 0);

guidata(hObject, handles);

function edit_winEdge_CreateFcn(hObject, eventdata, handles)

function push_pitchAnal_Callback(hObject, eventdata, handles)
%Find which syllable the popup is set to
ind = get(handles.popup_sylSelect,'Value');
renditions = length(handles.data.sylSnips{ind}.ts);
P.sr = handles.fs;
P.minf0 = 240;
P.maxf0 = 1200;
    
for i = 1:renditions
    %pnts = getWarpedStarts(handles.data.sylSnips{ind}.path{i},handles.data.sylSnips{ind}.winEdge);
    pnts = interp1(handles.data.sylSnips{ind}.path{i}(:,1),handles.data.sylSnips{ind}.path{i}(:,2),handles.data.sylSnips{ind}.winEdge);
    P.range = round((pnts*handles.fs/1000));

    R=yin(handles.data.sylSnips{ind}.ts{i},P);
    
    % Report f0 of the entire signal
	[~, idx] = min(R.ap0);
	best=R.f0(idx);
    handles.data.sylSnips{ind}.pitch(i) = (2^best)*440; %Convert octave into Hz
end

figure;
hist(handles.data.sylSnips{ind}.pitch,20);

guidata(hObject, handles);

function popup_sylSelect_Callback(hObject, eventdata, handles)

function popup_sylSelect_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_temporal1_Callback(hObject, eventdata, handles)

%Generate Workset given current selections
[handles.curAudio handles.curElements] = getWorkset(handles);

%Calculate Syllable and Gap Duration Distributions
[handles.data.temporal.sylDur,handles.data.temporal.gapDur,handles.data.temporal.sylDurbyType,handles.data.temporal.sylTypes] = SylGapDistr(handles.curElements);

%Plot duration distributions on the same axis
minval = min(min(handles.data.temporal.sylDur),min(handles.data.temporal.gapDur));
maxval = max(max(handles.data.temporal.sylDur),max(handles.data.temporal.gapDur));
bins = 100;
[sX,sY] = epdf(handles.data.temporal.sylDur,bins,minval,maxval);
[gX,gY] = epdf(handles.data.temporal.gapDur,bins,minval,maxval);

axes(handles.axes_SylGapDistr)
cla
hold on
plot(sX,sY,'b')
plot(gX,gY,'r')
axis tight
ylim([0,max([sY,gY])*1.1])
set(handles.axes_SylGapDistr,'XTick',0:0.1:(maxval+.05))
set(handles.axes_SylGapDistr,'YTick',0:0.1:max([sY,gY]))
hold off

%Plot syllable durations by type
axes(handles.axes_SylIntDur)
cla
hold on
for i = 1:length(handles.data.temporal.sylTypes)
   jitter = randn(1,length(handles.data.temporal.sylDurbyType{i}))*0.1;
   xs = (ones(1,length(handles.data.temporal.sylDurbyType{i})).*i)+jitter;
   scatter(xs,handles.data.temporal.sylDurbyType{i},'.')
   %scatter(handles.data.temporal.sylDurbyType{i},xs'.')
   tempTypes{i} = num2str(handles.data.temporal.sylTypes(i));
end

set(handles.axes_SylIntDur,'XTickLabel',tempTypes,'XTick',1:length(tempTypes))
xlim([0,length(tempTypes)+1])
set(handles.axes_SylIntDur,'XTickLabel',tempTypes,'XTick',1:length(tempTypes))
set(handles.axes_SylIntDur,'YTick',0:0.1:max(handles.data.temporal.sylDur))
hold off

guidata(hObject, handles);

function push_temporal2_Callback(hObject, eventdata, handles)
%Generate Workset given current selections
[handles.curAudio handles.curElements] = getWorkset(handles);

%Calculate beat/rhythm for the current workset
[handles.data.temporal.rhythm.mean,handles.data.temporal.rhythm.std] = CalcRhythm(handles.curAudio);

%Plot beat/rhythm spectrum
axes(handles.axes_rhythm)
cla
hold on
imagesc(fliplr(handles.data.temporal.rhythm.mean(5:40))')
axis tight; axis ij
set(handles.axes_rhythm,'XTick',[],'XTickLabel',{});
set(handles.axes_rhythm,'YTickLabel',{'40','30','20','10'},'YTick',1:10:40)
hold off

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Syntax Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_Syntax1_Callback(hObject, eventdata, handles)

%Generate Workset given current selections
[handles.curAudio handles.curElements] = getWorkset(handles);

%Generate sequence transition matrix
[handles.data.syntax.transMatrix,sylTypes] = SequenceMatrix(handles.curElements);

%Plot matrix
for i = 1:length(sylTypes)
    tempTypes{i} = num2str(sylTypes(i));
end
axes(handles.axes_transMatrix);
imagesc(handles.data.syntax.transMatrix);
axis tight
set(handles.axes_transMatrix,'XTickLabel',tempTypes,'XTick',1:length(sylTypes))
set(handles.axes_transMatrix,'YTickLabel',tempTypes,'YTick',1:length(sylTypes))
guidata(hObject, handles);

function push_Syntax2_Callback(hObject, eventdata, handles)

%Generate sequence transition matrix
[handles.data.syntax.transMatrix,sylTypes] = SequenceMatrixSE(handles.curElements);

h = figure;
axes1 = axes('Parent',h);%'YDir','reverse','XTick',[1 2 3 4 5],'Layer','top');
imagesc(handles.data.syntax.transMatrix);
hold on
axis tight

%Plot matrix
for i = 1:length(sylTypes)
    tempTypes{i} = num2str(sylTypes(i));
    totFol = sum(handles.data.syntax.transMatrix(i,:));
     for j = 1:length(sylTypes)
         text(j-.1,i,num2str(handles.data.syntax.transMatrix(i,j)/totFol,2));
     end
end
tempTypes{sylTypes == -99} = 'S';
tempTypes{sylTypes == 99} = 'E';
set(axes1,'XTickLabel',tempTypes,'XTick',1:length(sylTypes))
set(axes1,'YTickLabel',tempTypes,'YTick',1:length(sylTypes))
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_exportFeatures_Callback(hObject, eventdata, handles)

function push_exportTransMat_Callback(hObject, eventdata, handles)

%Generate sequence transition matrix
[handles.data.syntax.transMatrix,sylTypes] = SequenceMatrix(handles.curElements);

h = figure;
axes1 = axes('Parent',h);%'YDir','reverse','XTick',[1 2 3 4 5],'Layer','top');
imagesc(handles.data.syntax.transMatrix);
hold on
axis tight

%Plot matrix
for i = 1:length(sylTypes)
    tempTypes{i} = num2str(sylTypes(i));
    totFol = sum(handles.data.syntax.transMatrix(i,:));
    for j = 1:length(sylTypes)
        text(j-.1,i,num2str(handles.data.syntax.transMatrix(i,j)/totFol,2));
    end
end
set(axes1,'XTickLabel',tempTypes,'XTick',1:length(sylTypes))
set(axes1,'YTickLabel',tempTypes,'YTick',1:length(sylTypes))
hold off

function push_exportSongs_Callback(hObject, eventdata, handles)

function push_exportAnalysis_Callback(hObject, eventdata, handles)
