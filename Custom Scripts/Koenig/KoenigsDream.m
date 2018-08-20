function varargout = KoenigsDream(varargin)
% KOENIGSDREAM MATLAB code for KoenigsDream.fig
%      Last modified by TMO on 3/31/14
%      KOENIGSDREAM, by itself, creates a new KOENIGSDREAM or raises the existing
%      singleton*.
%
%      H = KOENIGSDREAM returns the handle to a new KOENIGSDREAM or the handle to
%      the existing singleton*.
%
%      KOENIGSDREAM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KOENIGSDREAM.M with the given input arguments.
%
%      KOENIGSDREAM('Property','Value',...) creates a new KOENIGSDREAM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before KoenigsDream_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to KoenigsDream_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help KoenigsDream

% Last Modified by GUIDE v2.5 20-Jun-2014 12:07:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @KoenigsDream_OpeningFcn, ...
                   'gui_OutputFcn',  @KoenigsDream_OutputFcn, ...
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

function KoenigsDream_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for KoenigsDream
handles.output = hObject;

%Initialize variables
handles.fs=44150;
handles.dirname =[];
handles.datasetFilename =[];
handles.elements = [];
handles.totFilelist = [];
handles.filelist =[];
handles.dataset =[];

%Clear single-shot axes
set(handles.axes_rendSpect,'XTick',[],'YTick',[]);
set(handles.axes_stackPlot,'XTick',[],'YTick',[]);
set(handles.axes_timeDistr,'XTick',[],'YTick',[]);
set(handles.axes_alignPaths,'XTick',[],'YTick',[]);
set(handles.axes_SylGapDistr,'XTick',[],'YTick',[]);
set(handles.axes_featStats,'XTick',[],'YTick',[]);
set(handles.axes_Freq,'XTick',[],'YTick',[]);
set(handles.axes_varDecomp,'XTick',[],'YTick',[]);
set(handles.axes_pitchSnips,'XTick',[],'YTick',[]);

%Clear longitudinal axes
set(handles.axes_longTemporal,'XTick',[],'YTick',[]);
set(handles.axes_longVar,'XTick',[],'YTick',[]);
set(handles.axes_longSpectral,'XTick',[],'YTick',[]);
set(handles.axes_longOptional,'XTick',[],'YTick',[]);


%Set initial values for feature popup boxes
handles.featureSet = {'Amp Mod';'Freq Mod';'Entropy';'1/2 Amp';'Grav Cent';...
    'Pitch Good'; 'Pitch'};
set(handles.popup_featStats,'String',handles.featureSet);

% Update handles structure
guidata(hObject, handles);

function varargout = KoenigsDream_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Independent, in-GUI functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles] = clearSingleShotAxes(handles)
% Clears out all of the axes that show changeable data
axes(handles.axes_rendSpect)
cla;
axes(handles.axes_stackPlot)
cla;
axes(handles.axes_timeDistr)
cla;
axes(handles.axes_alignPaths)
cla;
axes(handles.axes_SylGapDistr)
cla;
axes(handles.axes_featStats)
cla;
axes(handles.axes_pitchSnips)
cla;
axes(handles.axes_Freq)
cla;

function [handles] = clearLongitudinalAxes(handles)
% Clears out all of the axes that show changeable data
axes(handles.axes_longTemporal)
cla;
axes(handles.axes_longSpectral)
cla;
axes(handles.axes_longOptional)
cla;

function [handles] = SylGapDistr(handles,ind)
x = false;
if get(handles.check_durSylDistr,'Value')
    plot(handles.procSSData{ind}.sX,handles.procSSData{ind}.sY,'b')
    x = true;
end
if get(handles.check_durGapDistr,'Value')
    plot(handles.procSSData{ind}.gX,handles.procSSData{ind}.gY,'r')
    x = true;
end
if get(handles.check_durAllDistr,'Value')
    plot(handles.procSSData{ind}.tX,handles.procSSData{ind}.tY,'k')
    x = true;
end

y = false;
if get(handles.check_durBySyllGap,'Value')
    for i = 1:length(handles.dataset{ind}.sgLabels)
        jitter = randn(1,size(handles.procSSData{ind}.intervalMatrix,1))*0.1;
        xs = (ones(1,size(handles.procSSData{ind}.intervalMatrix,1)).*i)+jitter;
        scatter(xs,handles.procSSData{ind}.intervalMatrix(:,i),'.')
    end
    errorbar(1:length(handles.dataset{ind}.sgLabels),handles.procSSData{ind}.intervalMeans,handles.procSSData{ind}.intervalStd,'ok','LineWidth',1.5);
    set(handles.axes_SylGapDistr,'XTickLabel',handles.dataset{ind}.sgLabels,'XTick',1:length(handles.dataset{ind}.sgLabels))
    set(handles.axes_SylGapDistr,'YTick',0:50:max(handles.procSSData{ind}.intervalMatrix(:)))
    xlim([0,length(handles.dataset{ind}.sgLabels)+1])
    ylim([0,(max(handles.procSSData{ind}.intervalMatrix(:))+10)])
    y = true;
elseif get(handles.check_durBySyllNGap,'Value')
    for i = 1:length(handles.dataset{ind}.sNgLabels)
        jitter = randn(1,size(handles.procSSData{ind}.intervalSNGMatrix,1))*0.1;
        xs = (ones(1,size(handles.procSSData{ind}.intervalSNGMatrix,1)).*i)+jitter;
        scatter(xs,handles.procSSData{ind}.intervalSNGMatrix(:,i),'.')
    end
    errorbar(1:length(handles.dataset{ind}.sNgLabels),handles.procSSData{ind}.intervalSNGMeans,handles.procSSData{ind}.intervalSNGStd,'ok','LineWidth',1.5);
    set(handles.axes_SylGapDistr,'XTickLabel',handles.dataset{ind}.sNgLabels,'XTick',1:length(handles.dataset{ind}.sNgLabels))
    set(handles.axes_SylGapDistr,'YTick',0:50:max(handles.procSSData{ind}.intervalSNGMatrix(:)))
    xlim([0,length(handles.dataset{ind}.sNgLabels)+1])
    ylim([0,(max(handles.procSSData{ind}.intervalSNGMatrix(:))+10)])
    y = true;
end

if x
    xlim([0,max(handles.procSSData{ind}.tX)*1.1])
    ylim([0,max([handles.procSSData{ind}.tY,handles.procSSData{ind}.sY,handles.procSSData{ind}.gY])*1.1])
    set(handles.axes_SylGapDistr,'XTickLabel',0:50:(max(handles.procSSData{ind}.tX)+10),'XTick',0:50:(max(handles.procSSData{ind}.tX)+10))
    set(handles.axes_SylGapDistr,'YTick',0:0.1:max([handles.procSSData{ind}.tY,handles.procSSData{ind}.sY,handles.procSSData{ind}.gY]))
    set(handles.axes_SylGapDistr,'TickDir','out')
    xlabel('Interval Dur (ms)');
    ylabel('P(x)')
end

if y
    set(handles.axes_SylGapDistr,'TickDir','out')
    xlabel('Int Type');
    ylabel('Duration (ms)')
end

function [info] = getPitchSnipInfo(init, labels)
%Setup the dialog
numSyl = size(labels,1);
dlg_title = 'Input pitch snipping parameters';
num_lines = 1;
options.Resize='on';

prompt = []; def = [];
for i = 1:numSyl
    %Create labels for the dialog box
    line = {['Syllable ' char(labels(i)) ' offset'], ['Syllable ' char(labels(i)) ' center frequency']};
    prompt = [prompt, line];

    %Use either existing or default values
    if isempty(init)
        line = {'20','570'};
    else
        line = {num2str(init((2*i)-1)), num2str(init(2*i))};
    end
    def = [def, line];
end

%Present dialog
answer = inputdlg(prompt,dlg_title,num_lines,def,options);

%Parse dialog answers into numbers
for i = 1:size(answer,1)
    info(i) = str2num(answer{i});
end

function [filtAudio] = Prep(audio)
%Constants for Bandpass Audio (300-10000kHz)
HP_fNorm = 300/(44150/2);
LP_fNorm = 6500/(44150/2); %Changed 05/16/13 to match Farhan and Cengiz
[BP_b,BP_a] = butter(2,[HP_fNorm LP_fNorm]);

%Cycle through all each rendition and zero-phase filter
renditions = size(audio,2);
parfor i = 1:renditions
    if ~isempty(audio)
        filtAudio{i} = filtfilt(BP_b,BP_a,audio{i});
    end
end

function SylFeatureSet = parseSpectral(rendFeatures,path,tempBreaks)
%Get the list of all fields included in the features structurs
names = fieldnames(rendFeatures);

%Loop through each of the annotated syllables, parse the feature vectors,
%and file in to the appropriate structure
for i = 1:size(tempBreaks,1)
    %Get indices and pointers to pull and push the right data
    [rendBreaks] = getWarpedStarts(path,tempBreaks(i,:));
    ind = rendBreaks(1):rendBreaks(2);

    %Work through each field to parse the data
    for j = 1:length(names)
        %For ease of reading, copy rendition-length feature to tmp var
        tmp = eval(['rendFeatures.' names{j}]);
        eval(['SylFeatureSet.ts_' names{j} '{i} = tmp(ind);']);
        eval(['SylFeatureSet.mean_' names{j} '{i} = mean(tmp(ind));']);
    end
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

function [hourVect,minVect] = extractTimeVects(filenames)
%Extract the hour and minute value from the filenames
if iscellstr(filenames)
    rendNum = size(filenames,2);
else
    rendNum = size(filenames,1);
end
for i = 1:rendNum
    sp = regexp(filenames{i},'_','split');
    hourVect(i) = str2num(sp{6});
    minVect(i) = str2num(sp{7});
end

function IntMat = calcIntervals(data)
%Extract all interval durations from the paths and templatebreaks

%Get the template syllable breaks from the passed file
a = transpose(data.templatesyllBreaks);
a = a(:);

%Interval Durations from the paths and template files
rendNum = length(data.audio);
IntMat = [];
path = {};
warpedOut = {};
% for i = 1:rendNum
%     %Rover rendition path and rendition length
%     path = [data.p{i},data.q{i}];
%     
%     %Calculate intervals and add to the stack
%     [warpedOut] = getWarpedStarts(path,a(:));
%     IntMat = [IntMat; diff(warpedOut)'];
% end

parfor i = 1:rendNum
    %Rover rendition path and rendition length
    path{i} = [data.p{i},data.q{i}];
    
    %Calculate intervals and add to the stack
    [warpedOut{i}] = getWarpedStarts(path{i},a(:));
    IntMat(i,:) = diff(warpedOut{i})';
end

function h = updateTimeDistr(handles,hourVect,minVect)
%Break down the matched hour and minute vectors to display a histogram of
%the recording times. Axes are specified in the parent function.

%Width of bins in hours; you must change the 'bar' x-limits to match,
%otherwise it'll throw an error.
binWidth = 0.25; 
ToD = hourVect+minVect/60;
hs = histc(ToD,0:binWidth:24);

h = bar(.25:binWidth:23.75,hs(2:end-1));

function h = updatePathsPlot(p,q)
%Loop through the passed alignment path cell arrays and plot
hold on
for i = 1:length(p)
    cumPath(i,:) = q{i}-p{i};
    h = plot(p{i},cumPath(i,:));
end
h = plot(mean(cumPath,1),'k','LineWidth',2);
hold off

function [featMat,featMeans,featStd] = updatefeatStats(syl_features,feat,labels)

%Generate feature arrays for selected feature
for i = 1:length(syl_features)
    featMat(i,:) = cell2mat(eval(['syl_features{i}.' char(feat)]));
end

%Calculate descriptive stats
featMeans = mean(featMat,1);
featStd = std(featMat,1);

%Scatter plot results
for i = 1:length(labels)
    jitter = randn(1,size(featMat,1))*0.1;
    xs = (ones(1,size(featMat,1)).*i)+jitter;
    scatter(xs,featMat(:,i),'.')
end
errorbar(1:length(labels),featMeans,featStd,'ok','LineWidth',1.5)

function updatePitchSnips(pitchSnips, labels)
%Calculate descriptive stats
pitchMeans = mean(pitchSnips,1);
pitchStd = std(pitchSnips,1);

%Scatter plot results
for i = 1:length(labels)
    jitter = randn(1,size(pitchSnips,1))*0.1;
    xs = (ones(1,size(pitchSnips,1)).*i)+jitter;
    scatter(xs,pitchSnips(:,i),'.')
end
errorbar(1:length(labels), pitchMeans, pitchStd,'ok','LineWidth',1.5)

function handles = updateSimStats(handles)
%This is the number of syllables to be randomly chosen for the "quick"
%analysis
quickNum = 20;

%The amount of jitter/slack (in seconds) allowable in MatchScore
slack = 0.01;

%Generate feature arrays for each feature
axes(handles.axes_Freq)
set(handles.axes_Freq,'NextPlot','ReplaceChildren')
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

set(handles.axes_Freq,'FontSize',10)
xlim([0,length(handles.data.spectral.syl_featuresKEY)+1])
ylim([-.1,1.1])
set(handles.axes_Freq,'XTickLabel',tempTypes,'XTick',1:(length(handles.data.spectral.syl_featuresKEY)+1))
set(handles.axes_Freq,'YTickLabel',ylabel,'YTick',ytick)

hold off

function [procSSData, pitchSnipParams] = processDataset(handles,ind)
%This function manages all data processing for the GUI. It's used in both
%the single shot analysis and the longitudinal analysis to minimize
%coding duplication.

%Copy over the dataset name to maintain chain back to data
procSSData.datasetName = handles.datatitles{ind};

%Capture the time of recording from the filename (all of Tim's recording
%programs are consistent with this format.)
[procSSData.hourVect,procSSData.minVect] = extractTimeVects(handles.dataset{ind}.renditionFilenames);

%Calculate interval lengths for syls and gaps together and apart
procSSData.intervalMatrix = calcIntervals(handles.dataset{ind});
m = 1;
for i = 1:2:(size(procSSData.intervalMatrix,2)-1)
    procSSData.intervalSNGMatrix(:,m) = procSSData.intervalMatrix(:,i)+procSSData.intervalMatrix(:,i+1);
    m = m+1;
end
procSSData.intervalSNGMatrix(:,m) = procSSData.intervalMatrix(:,end);

%Calculate the interval duration distributions
minval = min(procSSData.intervalMatrix(:))-10;
maxval = max(procSSData.intervalMatrix(:))+10;
bins = 100;
[procSSData.sX,procSSData.sY] = epdf(procSSData.intervalMatrix(:,1:2:end),bins,minval,maxval);
[procSSData.gX,procSSData.gY] = epdf(procSSData.intervalMatrix(:,2:2:end),bins,minval,maxval);
[procSSData.tX,procSSData.tY] = epdf(procSSData.intervalMatrix,bins,minval,maxval);

%Calculate descriptive statistic
procSSData.intervalMeans = mean(procSSData.intervalMatrix,1);
procSSData.intervalStd = std(procSSData.intervalMatrix,1);

procSSData.intervalSNGMeans = mean(procSSData.intervalSNGMatrix,1);
procSSData.intervalSNGStd = std(procSSData.intervalSNGMatrix,1);

%Calculate the features for the dataset
path = {};
parfor i = 1:length(handles.dataset{ind}.audio)
    %Get spectral features for the entire motif
    rend_features{i} = koenigSpectral(handles.dataset{ind}.audio{i},handles.fs);

    %Parse the feature arrays by syllable boundaries
    path{i} = [handles.dataset{ind}.p{i},handles.dataset{ind}.q{i}];
    syl_features{i} = parseSpectral(rend_features{i},path{i},handles.dataset{ind}.templatesyllBreaks);

    %Calculate spectral variability byt some yet to be named measure

end
procSSData.rend_features = rend_features;
procSSData.syl_features = syl_features;

%Calculate the Pitch Snip for each syllable (given the user specified offset)
if ~isfield(handles,'pitchSnipParams')
    handles.pitchSnipParams = getPitchSnipInfo([], handles.dataset{ind}.sLabels);
end
procSSData.pitchSnips = pitchSnip(handles.dataset{ind}.audio, handles.dataset{ind}.p, handles.dataset{ind}.q, handles.dataset{ind}.templatesyllBreaks, handles.pitchSnipParams);
pitchSnipParams = handles.pitchSnipParams;


function pitchCell = pitchTraces(audio, p, q,  templatesyllBreaks, target, centerFreq, buffer)
tic
%Select target syllable bounds
a = templatesyllBreaks(target,:);

%Pitch Calc Parameters
fs = 44150; %Recording sampling rate
winSize = 220; %Length (in samples) of syllable segment to analyze (5ms)
winStep = 44; %window advance (1ms)
numHarmonics = 20; %Number of harmonics to locate

%Interval Durations from the paths and template files
rendNum = length(audio);
pitchCell = {};
warpedOut = [];
for i = 1:rendNum
    %Rover rendition path and rendition length
    path{i} = [p{i},q{i}];

    %Find the target interval position
    warpedOut = getWarpedStarts(path{i},a);

    %Get audio snippet
    startPnt(i) = floor(fs*(warpedOut(1)/1000)) - buffer;
    endPnt(i)  = ceil(fs*(warpedOut(2)/1000)) + buffer;
    
    %Window the snippet
    windowedSnip{i} = window(audio{i}(startPnt(i):endPnt(i))', winSize, winStep);

    %Cycle through the snippet windows and calculate pitch
    pitch = [];
    parfor j = 1:size(windowedSnip{i},1)
        %Calculate pitch
        [pitch(j), ~] = sparse_fftm(windowedSnip{i}(j,:), centerFreq, numHarmonics, fs);
    end
    pitchCell = [pitchCell; pitch];
end
toc


function pitchMat = pitchSnip(audio, p, q,  templatesyllBreaks, options)

%Select target syllable bounds
a = transpose(templatesyllBreaks(:,1));

%Pitch Calc Parameters
pLength = 220; %Length (in samples) of syllable segment to analyze
numHarmonics = 20; %Number of harmonics to locate
fs = 44150; %Recording sampling rate

%Interval Durations from the paths and template files
rendNum = length(audio);
numSyls = size(templatesyllBreaks,1);
pitchMat = [];
for i = 1:rendNum
    %Rover rendition path and rendition length
    path{i} = [p{i},q{i}];
    
    %Calculate pitch and add to the stack
    [warpedOut(i,:)] = getWarpedStarts(path{i},a(:));
    
    %Step through each syllable and calc pitch, or give NaN
    parfor j = 1:numSyls
            %Get syllable targeting info from options file
            pOffset(j) = round((options((2*j)-1)/1000)*fs);
            centerFreq(j) = options(2*j);
            
            %Get audio snippet
            startPnt(j) = round(fs*(warpedOut(i,j)/1000))+pOffset(j);
            endPnt(j)  = min(length(audio{i}),(startPnt(j)+pLength-1));
            snip(j,:) = audio{i}(startPnt(j):endPnt(j));
            
            %Calculate pitch
            [pitchRend(j), ~] = sparse_fftm(snip(j,:),centerFreq(j),numHarmonics,fs);
    end
    pitchMat = [pitchMat; pitchRend];
end

function handles = updateLongTemporal(handles)
%This function is called as the timeseries selections for the longitudinal temporal plots are
%changed. It's also called during the export image request.

%Clear all arrays to prevent confusion
handles.long.temporalBlue = [];
handles.long.temporalRed = [];
handles.long.temporalGreen = [];
xmin = inf;
xmax = -inf;
ymin = inf;
ymax = -inf;

%Collect syl/gap type and interval selections
intType = get(handles.check_longSylGap,'Value');
blueType = get(handles.popup_temporalBlue,'Value');
redType = get(handles.popup_temporalRed,'Value');
greenType = get(handles.popup_temporalGreen,'Value');
contents = get(handles.popup_temporalBlue,'String');

%Extract values for the blue trace (if it's selected)
handles.long.temporalBlue.type = contents(blueType);
if ~strcmp(handles.long.temporalBlue.type, 'NA')
    if strcmp(handles.long.temporalBlue.type, 'Total')
        [handles.long.temporalBlue.time,handles.long.temporalBlue.mean,handles.long.temporalBlue.std,handles.long.temporalBlue.name] = extractLongTime(handles.procSSData,-1,blueType);
    else
         [handles.long.temporalBlue.time,handles.long.temporalBlue.mean,handles.long.temporalBlue.std,handles.long.temporalBlue.name] = extractLongTime(handles.procSSData,intType,blueType);
    end
    
    xmin = min(xmin,min(handles.long.temporalBlue.time)); xmax = max(xmax,max(handles.long.temporalBlue.time));
    ymin = min(ymin,min(handles.long.temporalBlue.mean-handles.long.temporalBlue.std)); ymax = max(ymax,max(handles.long.temporalBlue.mean+handles.long.temporalBlue.std));
end

%Extract values for the red trace (if it's selected)
handles.long.temporalRed.type = contents(redType);
if ~strcmp(handles.long.temporalRed.type, 'NA')
    if strcmp(handles.long.temporalRed.type, 'Total')
        [handles.long.temporalRed.time,handles.long.temporalRed.mean,handles.long.temporalRed.std,handles.long.temporalRed.name] = extractLongTime(handles.procSSData,-1,redType);
    else
         [handles.long.temporalRed.time,handles.long.temporalRed.mean,handles.long.temporalRed.std,handles.long.temporalRed.name] = extractLongTime(handles.procSSData,intType,redType);
    end
    
    xmin = min(xmin,min(handles.long.temporalRed.time)); xmax = max(xmax,max(handles.long.temporalRed.time));
    ymin = min(ymin,min(handles.long.temporalRed.mean-handles.long.temporalRed.std)); ymax = max(ymax,max(handles.long.temporalRed.mean+handles.long.temporalRed.std));

end

%Extract values for the green trace (if it's selected)
handles.long.temporalGreen.type = contents(greenType);
if ~strcmp(handles.long.temporalGreen.type, 'NA')
    if strcmp(handles.long.temporalGreen.type, 'Total')
        [handles.long.temporalGreen.time,handles.long.temporalGreen.mean,handles.long.temporalGreen.std,handles.long.temporalGreen.name] = extractLongTime(handles.procSSData,-1,greenType);
    else
         [handles.long.temporalGreen.time,handles.long.temporalGreen.mean,handles.long.temporalGreen.std,handles.long.temporalGreen.name] = extractLongTime(handles.procSSData,intType,greenType);
    end
    
    xmin = min(xmin,min(handles.long.temporalGreen.time)); xmax = max(xmax,max(handles.long.temporalGreen.time));
    ymin = min(ymin,min(handles.long.temporalGreen.mean-handles.long.temporalGreen.std)); ymax = max(ymax,max(handles.long.temporalGreen.mean+handles.long.temporalGreen.std));

end

%Plot all of the selected functions
axes(handles.axes_longTemporal)
cla; hold on
if ~strcmp(handles.long.temporalBlue.type, 'NA')
    if ~get(handles.check_tempNormal,'Value')
        plot(handles.long.temporalBlue.time,handles.long.temporalBlue.mean,'b')
        errorbar(handles.long.temporalBlue.time,handles.long.temporalBlue.mean,handles.long.temporalBlue.std,'xb')
    else
        plot(handles.long.temporalBlue.time,handles.long.temporalBlue.mean/handles.long.temporalBlue.mean(1),'b')
        errorbar(handles.long.temporalBlue.time,handles.long.temporalBlue.mean/handles.long.temporalBlue.mean(1),handles.long.temporalBlue.std/handles.long.temporalBlue.mean(1),'xb')
    end
    
%     if get(handles.check_tempCV,'Value')
%         plot(handles.long.temporalBlue.time,handles.long.temporalBlue.std./handles.long.temporalBlue.mean,'ob')
%     end
end
if ~strcmp(handles.long.temporalRed.type, 'NA')
    if ~get(handles.check_tempNormal,'Value')
        plot(handles.long.temporalRed.time,handles.long.temporalRed.mean,'r')
        errorbar(handles.long.temporalRed.time,handles.long.temporalRed.mean,handles.long.temporalRed.std,'xr')
    else
        plot(handles.long.temporalRed.time,handles.long.temporalRed.mean/handles.long.temporalRed.mean(1),'r')
        errorbar(handles.long.temporalRed.time,handles.long.temporalRed.mean/handles.long.temporalRed.mean(1),handles.long.temporalRed.std/handles.long.temporalRed.mean(1),'xr')        
    end
end
if ~strcmp(handles.long.temporalGreen.type, 'NA')
    if ~get(handles.check_tempNormal,'Value')
        plot(handles.long.temporalGreen.time,handles.long.temporalGreen.mean,'g')
        errorbar(handles.long.temporalGreen.time,handles.long.temporalGreen.mean,handles.long.temporalGreen.std,'xg')
    else
        plot(handles.long.temporalGreen.time,handles.long.temporalGreen.mean/handles.long.temporalGreen.mean(1),'g')
        errorbar(handles.long.temporalGreen.time,handles.long.temporalGreen.mean/handles.long.temporalGreen.mean(1),handles.long.temporalGreen.std/handles.long.temporalGreen.mean(1),'xg')        
    end
end
hold off; box off

%Format the axes
xlim([(xmin-0.2*(xmax-xmin)),(xmax+0.2*(xmax-xmin))])
xrangeVect = round(xmin:xmax);
if length(xrangeVect) > 6
    xT = round(xrangeVect(round(linspace(1,length(xrangeVect),5))));
else
    xT = xrangeVect;
end



ylim([ymin, ymax])
yrangeVect = round(ymin:ymax);
yT = round(linspace(yrangeVect(1),yrangeVect(end),5)*100)/100;
if get(handles.check_tempNormal,'Value')
    ylim([0.5,1.5])
    yT = 0.5:0.25:1.5;
end

set(handles.axes_longTemporal,'XTickLabel',xT,'XTick',xT)
set(handles.axes_longTemporal,'YTickLabel',yT,'YTick',yT)
set(handles.axes_longTemporal,'TickDir','out')
xlabel('Time (days)');
ylabel('Duration (ms)')

function [timepnts,meanVal,stdVal,filename] = extractLongTime(procSSData,intType,intNum)
%From each array within procSSData, this function extracts the duration information for the interval
%selected by the user. Data is returned as structured arrays to maximize plotting efficiency.
%Initialize variables
timepnts = [];
meanVal = [];
stdVal = [];

%Cycle through each array in procSSData to extract the required values
numDatasets = length(procSSData);

for i = 1:numDatasets
    if intType == 1                 %Syls and Gaps
        meanVal(i) = procSSData{i}.intervalMeans (intNum-1);
        stdVal(i) = procSSData{i}.intervalStd(intNum-1);
    elseif intType == 0         %Syls + Gaps
        meanVal(i) = procSSData{i}.intervalSNGMeans(intNum-1);
        stdVal(i) = procSSData{i}.intervalSNGStd(intNum-1);
    elseif intType == -1        %Total Length
        meanVal(i) = mean(sum(procSSData{i}.intervalMatrix,2));
        stdVal(i) = std(sum(procSSData{i}.intervalMatrix,2));
    end
    timepnts (i) = procSSData{i}.recDateNum;
    filename{i} = procSSData{i}.datasetName;
end

%Sort data by timepoints, but keep everything aligned
[timepnts, IX] = sort(timepnts);
meanVal = meanVal(IX);
stdVal = stdVal(IX);
filename = filename{IX};

function handles = updateLongSpectral(handles)
%This function is called as the timeseries selections for the longitudinal spectral plots are
%changed. It's also called during the export image request.

%Clear all arrays to prevent confusion
handles.long.spectralBlue = [];
handles.long.spectralRed = [];
handles.long.spectralGreen = [];
xmin = inf; xmax = -inf;
ymin = inf; ymax = -inf;

%Collect syl/gap type and interval selections
featType = get(handles.popup_featSelect,'Value');
blueType = get(handles.popup_spectralBlue,'Value');
redType = get(handles.popup_spectralRed,'Value');
greenType = get(handles.popup_spectralGreen,'Value');
contents = get(handles.popup_spectralBlue,'String');
featContents = get(handles.popup_featSelect,'String');

%Extract values for the blue trace (if it's selected)
handles.long.spectralBlue.type = contents(blueType);
handles.long.spectralBlue.feat = featContents(featType);
if ~strcmp(handles.long.spectralBlue.type, 'NA')
     option = [];
    [handles.long.spectralBlue.time,handles.long.spectralBlue.mean,handles.long.spectralBlue.std,handles.long.spectralBlue.name] = extractLongSpectral(handles.procSSData,featType,blueType,option);
    
    xmin = min(xmin,min(handles.long.spectralBlue.time)); xmax = max(xmax,max(handles.long.spectralBlue.time));
    ymin = min(ymin,min(handles.long.spectralBlue.mean-handles.long.spectralBlue.std)); ymax = max(ymax,max(handles.long.spectralBlue.mean+handles.long.spectralBlue.std));
end

%Extract values for the Red trace (if it's selected)
handles.long.spectralRed.type = contents(redType);
handles.long.spectralRed.feat = featContents(featType);
if ~strcmp(handles.long.spectralRed.type, 'NA')
     option = [];
    [handles.long.spectralRed.time,handles.long.spectralRed.mean,handles.long.spectralRed.std,handles.long.spectralRed.name] = extractLongSpectral(handles.procSSData,featType,redType,option);
    
    xmin = min(xmin,min(handles.long.spectralRed.time)); xmax = max(xmax,max(handles.long.spectralRed.time));
    ymin = min(ymin,min(handles.long.spectralRed.mean-handles.long.spectralRed.std)); ymax = max(ymax,max(handles.long.spectralRed.mean+handles.long.spectralRed.std));
end

%Extract values for the blue trace (if it's selected)
handles.long.spectralGreen.type = contents(greenType);
handles.long.spectralGreen.feat = featContents(featType);
if ~strcmp(handles.long.spectralGreen.type, 'NA')
     option = [];
    [handles.long.spectralGreen.time,handles.long.spectralGreen.mean,handles.long.spectralGreen.std,handles.long.spectralGreen.name] = extractLongSpectral(handles.procSSData,featType,greenType,option);
    
    xmin = min(xmin,min(handles.long.spectralGreen.time)); xmax = max(xmax,max(handles.long.spectralGreen.time));
    ymin = min(ymin,min(handles.long.spectralGreen.mean-handles.long.spectralGreen.std)); ymax = max(ymax,max(handles.long.spectralGreen.mean+handles.long.spectralGreen.std));
end

%Plot all of the selected functions
axes(handles.axes_longSpectral)
cla; hold on
if ~strcmp(handles.long.spectralBlue.type, 'NA')
    if ~get(handles.check_specNormal,'Value')
        plot(handles.long.spectralBlue.time,handles.long.spectralBlue.mean,'b')
        errorbar(handles.long.spectralBlue.time,handles.long.spectralBlue.mean,handles.long.spectralBlue.std,'xb')
    else
        plot(handles.long.spectralBlue.time,handles.long.spectralBlue.mean/handles.long.spectralBlue.mean(1),'b')
        errorbar(handles.long.spectralBlue.time,handles.long.spectralBlue.mean/handles.long.spectralBlue.mean(1),handles.long.spectralBlue.std/handles.long.spectralBlue.mean(1),'xb')
    end
end

if ~strcmp(handles.long.spectralRed.type, 'NA')
    if ~get(handles.check_specNormal,'Value')
        plot(handles.long.spectralRed.time,handles.long.spectralRed.mean,'r')
        errorbar(handles.long.spectralRed.time,handles.long.spectralRed.mean,handles.long.spectralRed.std,'xr')
    else
        plot(handles.long.spectralRed.time,handles.long.spectralRed.mean/handles.long.spectralRed.mean(1),'r')
        errorbar(handles.long.spectralRed.time,handles.long.spectralRed.mean/handles.long.spectralRed.mean(1),handles.long.spectralRed.std/handles.long.spectralRed.mean(1),'xr')
    end
end

if ~strcmp(handles.long.spectralGreen.type, 'NA')
    if ~get(handles.check_specNormal,'Value')
        plot(handles.long.spectralGreen.time,handles.long.spectralGreen.mean,'g')
        errorbar(handles.long.spectralGreen.time,handles.long.spectralGreen.mean,handles.long.spectralGreen.std,'xg')
    else
        plot(handles.long.spectralGreen.time,handles.long.spectralGreen.mean/handles.long.spectralGreen.mean(1),'g')
        errorbar(handles.long.spectralGreen.time,handles.long.spectralGreen.mean/handles.long.spectralGreen.mean(1),handles.long.spectralGreen.std/handles.long.spectralGreen.mean(1),'xg')
    end
end
hold off

%Format the axes
xlim([(xmin-0.2*(xmax-xmin)),(xmax+0.2*(xmax-xmin))])
xrangeVect = round(xmin:xmax);
xT = round(xrangeVect(round(linspace(1,length(xrangeVect),5))));

ylim([ymin, ymax])
yrangeVect = (ymin:ymax);
yT = (linspace(ymin,ymax,5));
if get(handles.check_specNormal,'Value')
    ylim([0.5,1.5])
    yT = 0.5:0.25:1.5;
end

set(handles.axes_longSpectral,'XTickLabel',xT,'XTick',xT)
set(handles.axes_longSpectral,'YTickLabel',yT,'YTick',yT)
set(handles.axes_longSpectral,'TickDir','out')
xlabel('Time (days)');
ylabel('Feature')

function [timepnts,meanVal,stdVal,filename] = extractLongSpectral(procSSData,feat,intNum,option)
%From each array within procSSData, this function extracts the selected information for the syllable
%selected by the user. Data is returned as structured arrays to maximize plotting efficiency.
%Initialize variables
timepnts = [];
meanVal = [];
stdVal = [];

%Cycle through each array in procSSData to extract the required values
numDatasets = length(procSSData);

for i = 1:numDatasets
    if feat  > 2                 %SAM features
        %Figure out what the user has asked to display
        featNames =  fieldnames(procSSData{i}.syl_features{1});
        k = strfind(featNames,'mean_');
        for j = 1:length(k)
            bools(j) = ~isempty(k{j});
        end
        meanFeats = featNames(bools);
        target = meanFeats(feat-2);
      
        %Generate feature arrays for selected feature
        for j = 1:length(procSSData{i}.syl_features)
            featMat(j,:) = cell2mat(eval(['procSSData{i}.syl_features{j}.' char(target)]));
        end

        %Calculate descriptive stats
        featMeans = mean(featMat,1);
        featStd = std(featMat,1);
        
        meanVal(i) = featMeans(intNum-1);
        stdVal(i) = featStd(intNum-1);
    elseif feat == 1         %Pitch Snips
         meanVal(i) = mean(procSSData{i}.pitchSnips(:,intNum-1));
         stdVal(i) = std(procSSData{i}.pitchSnips(:,intNum-1));
    elseif feat == 2        %SAP Var
%        meanVal(i) = mean(sum(procSSData{i}.intervalMatrix,2));
%        stdVal(i) = std(sum(procSSData{i}.intervalMatrix,2));
    end
    timepnts (i) = procSSData{i}.recDateNum;
    filename{i} = procSSData{i}.datasetName;
end

%Sort data by timepoints, but keep everything aligned
[timepnts, IX] = sort(timepnts);
meanVal = meanVal(IX);
stdVal = stdVal(IX);
filename = filename{IX};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function push_addDataset_Callback(hObject, eventdata, handles)
%This function adds new dataset files to the filelist, sorts them, and
%displays in the listbox. Some basic processing is done, but no data is
%loaded from the source. (That requires an additional button click.)

%Get the locations of the processed files (multiselect is possible)
[temp,path] = uigetfile('*.mat','Select the dataset file.','MultiSelect','on');

if isequal(temp,0) || isequal(path,0)
    %If they hit cancel of choose something dumb, generate error and end.
    set(handles.text_message,'String','Dataset location invalid or cancelled. Pick a valid file.');
else
    %If they don't, then parse the locations.
    fnames = [];
    if iscell(temp)
        %If there is more than one files selected, struct it
         for i = 1:size(temp,2)
             fnames{i} = char(temp{i});
             pathTot{i} = [path char(temp{i})];
         end
    else
        %Otherwise, just copy it out.
        fnames = temp;
        pathTot = [path,temp];
        i = 1;
    end

     if ~isempty(handles.datasetFilename)
         button = questdlg('Do you want to add this dataset(s) to the current list or replace them?','Add an annotation','Add to current','Replace current','Add to current');
         if (strcmp(button,'Replace current'))
            %Add here whatever other shit needs to be cleared from memory
            if isfield(handles,'dataset')
                handles = rmfield(handles,'dataset');
            end
            if isfield(handles,'procSSData')
                handles = rmfield(handles,'procSSData');
            end
            if isfield(handles, 'procCount')
                handles = rmfield(handles, 'procCount');
            end
            if isfield(handles,'long')
                handles = rmfield(handles,'long');
            end
            
            %Copy over the new files
             handles.datasetFilename = pathTot;
             handles.datatitles = fnames;

             %Sort the datatitles to alphabetically (which is also temporally if 
             %only one bird) and use the index to sort the datasetFilename as well
             if iscell(handles.datatitles)
                [handles.datatitles, ind] = sort(handles.datatitles);
                handles.datasetFilename = handles.datasetFilename(ind);
             end

             %Set the message to show annotation replaced
             set(handles.listbox_datasets,'String',handles.datatitles)
             set(handles.text_message,'String',['Added ' num2str(i) ' new datasets to the active list']);
             
         %This is the code for adding to the existing list 
         elseif (strcmp(button,'Add to current'))
             %Check if there are multiple files in queue already
             if iscell(handles.datasetFilename)
                 dFtemp = handles.datasetFilename{:};
                 dTtemp = handles.datatitles{:};
             else
                 dFtemp = handles.datasetFilename;
                 dTtemp = handles.datatitles;
             end
             
             %Check if there are multiple files to add
             if iscell(fnames)
                handles.datasetFilename = {dFtemp pathTot{:}}';
                handles.datatitles = {dTtemp fnames{:}}';
             else
                handles.datasetFilename = {dFtemp pathTot}';
                handles.datatitles = {dTtemp fnames}';                 
             end
             
             %Remove duplicates
             handles.datatitles = unique(handles.datatitles);
             handles.datasetFilename = unique(handles.datasetFilename);
                         
             %Sort the datatitles to alphabetically (which is also temporally if 
             %only one bird) and use the index to sort the datasetFilename as well
             if iscell(handles.datatitles)
                [handles.datatitles, ind] = sort(handles.datatitles);
                handles.datasetFilename = handles.datasetFilename(ind);
             end

             %Set the message to show annotation replaced
             set(handles.listbox_datasets,'String',handles.datatitles)
             set(handles.text_message,'String',['Added ' num2str(i) ' new datasets to the active list']);
         end
     %This is the code for loading the first/only folder    
     else
         handles.datasetFilename = pathTot;
         handles.datatitles = fnames;

         %Sort the datatitles alphabetically (which is also chronologically if 
         %only one bird) and use the index to sort the datasetFilename as well
         if iscell(handles.datatitles)
            [handles.datatitles, ind] = sort(handles.datatitles);
            handles.datasetFilename = handles.datasetFilename(ind);
         end

         %Set the message to show annotation replaced
         set(handles.listbox_datasets,'String',handles.datatitles)
         set(handles.text_message,'String',['Added ' num2str(i) ' new datasets to the active list']);
     end
     
end
guidata(hObject, handles);

function push_clearDatasets_Callback(hObject, eventdata, handles)
%Get user confirmation on clearing the list
button = questdlg('Are you sure you want to clear all datasets?','Clear datasets?','Clear''em All!','Nooooo!','Nooooo!');

if (strcmp(button,'Clear''em All!'))
    %Clear pointers to the datasets and their labels
     handles.datasetFilename = [];
     handles.datatitles = [];
     
    %Add here whatever other shit needs to be cleared from memory
    if isfield(handles,'dataset')
        handles = rmfield(handles,'dataset');
    end
    if isfield(handles,'procSSData')
        handles = rmfield(handles,'procSSData');
    end
    if isfield(handles, 'procCount')
        handles = rmfield(handles, 'procCount');
    end
    if isfield(handles,'long')
        handles = rmfield(handles,'long');
    end

    %Set the message to show datasets replaced
    set(handles.listbox_datasets,'String',handles.datatitles)
    set(handles.text_message,'String','All datasets cleared.');
else
    set(handles.text_message,'String','Datasets are unchanged.');
end

guidata(hObject, handles);

function push_loadDataset_Callback(hObject, eventdata, handles)
%Load the data from the saved .dat and .wav files to the workspace
if isempty(handles.datasetFilename)
    set(handles.text_message,'String','No datasets found to load. Check list, add files and try again.')
else
    %Add here whatever other shit needs to be cleared from memory
    if isfield(handles,'dataset')
        handles = rmfield(handles,'dataset');
    end
    if isfield(handles,'procSSData')
        handles = rmfield(handles,'procSSData');
    end
    if isfield(handles, 'procCount')
        handles = rmfield(handles, 'procCount');
    end
    if isfield(handles,'long')
        handles = rmfield(handles,'long');
    end
     
     %Work through each file to extract the audio data
    if iscell(handles.datasetFilename)
        numDatasets = size(handles.datasetFilename,2);
    else
        numDatasets = size(handles.datasetFilename,1);
    end
    
    %Grab the range of data to be pulled in for each file
    if get(handles.check_loadFull, 'Value')
        handles.parseRange = [0,24];
    else
        handles.parseRange = str2num(get(handles.edit_parseTimes,'String'));
    end
    [numParsed, col] = size(handles.parseRange);
    if col ~=2
        set(handles.text_message, 'String', 'Parse range incorrectly formatted. Correct and try again');
        return
    end
    handles.datatitles = {};
    setInx = 1;
    for i=1:numDatasets
        %Load specific variables from file
        if iscell(handles.datasetFilename)
            load(handles.datasetFilename{i},'data','filenames','sequence')
        else
            load(handles.datasetFilename,'data','filenames','sequence')
        end
        
        %Parse the filenames for the time of recordings
         [hourVect,minVect] = extractTimeVects(filenames);
         timeVect = hourVect + minVect./60;
        
        %Parse the sequence data to strip out the rendition count
        sp = regexp(sequence,' ','split');
        seqCheck(i) = str2double(sp(1));
        seq = str2double(sp(1));
        
        %Create Interval Labels
        sgLabels = {};
        sNgLabels = {};
        sLabels = {};
        for j = 1:length(sp{1})
            %Labels for breaking out gaps from syllables
            sgLabels = [sgLabels; sp{1}(j)];
            if j ~= length(sp{1}) 
                sgLabels = [sgLabels; [sp{1}(j) 'G']];
            end
            
            %Labels for joinging syllables and gaps
            if j ~= length(sp{1})
                sNgLabels = [sNgLabels; [sp{1}(j) '+G']];
            else
                sNgLabels = [sNgLabels; sp{1}(j)];
            end
            
            %Labels for just syllables
            sLabels = [sLabels; sp{1}(j)];
        end        
        
        %Parse the dataset date
        sp = regexp(filenames{1},'_','split');
        dateNum = datenum([sp{3} '-' sp{4} '-' sp{5}], 'yyyy-mm-dd');

        
        for k = 1:numParsed
            %Create the parsing index for this subset
            numMotif = str2double(get(handles.edit_numMotifs, 'String'));
            Indx = find(timeVect >= handles.parseRange(k,1) & timeVect <= handles.parseRange(k,2));
            if numMotif == -1
                rendIndx = Indx;
            else
                if length(Indx) > numMotif
                    rendIndx = Indx(1:numMotif);
                else
                    rendIndx = Indx;
                end
            end
            
            if ~isempty(rendIndx)
                %Create the name of the dataset snippet
                n = char(handles.datasetFilename(i));
                n0 = [n(1:end-4) ' ' num2str(handles.parseRange(k,1)) 'h'];
                handles.dataset{setInx}.datasetFilename = n0;
                n1 = regexp(n0, filesep, 'split');
                handles.datatitles = [handles.datatitles; n1(end)];


                handles.dataset{setInx}.renditionFilenames = filenames(rendIndx);

                %Parse the dataset date
                meanTime = mean([handles.parseRange(k,1), handles.parseRange(k,2)]);
                handles.dataset{setInx}.datenum = dateNum + (meanTime/24);

                %Parse the sequence data to strip out the rendition count
                handles.dataset{setInx}.sequence = seq;

                %Create Interval Labels
                handles.dataset{setInx}.sgLabels = sgLabels;
                handles.dataset{setInx}.sNgLabels = sNgLabels;
                handles.dataset{setInx}.sLabels = sLabels;

                if ~isempty(rendIndx)
                    %Copy over the actual audio data from alignment
                    handles.dataset{setInx}.audio = Prep(data.audio(rendIndx));
                    if isfield(data,'aligned_audioCube')
                        %handles.dataset{i}.audioCube = data.aligned_audioCube; %This should be changed so that I only keep the summed values for display in the stack plots... it's just wasteful
                        handles.dataset{setInx}.audioCube = -squeeze(sum(data.aligned_audioCube(rendIndx,:,:),2));
                    else
                        handles.dataset{setInx}.audioCube = [];
                    end
                    handles.dataset{setInx}.templatesyllBreaks = data.templatesyllBreaks;
                    handles.dataset{setInx}.p = data.p(rendIndx);
                    handles.dataset{setInx}.q = data.q(rendIndx);
                else
                    %Copy over the actual audio data from alignment
                    handles.dataset{setInx}.audio = {};
                    if isfield(data,'aligned_audioCube')
                        handles.dataset{setInx}.audioCube = {};
                    else
                        handles.dataset{setInx}.audioCube = {};
                    end
                    handles.dataset{setInx}.templatesyllBreaks = data.templatesyllBreaks;
                    handles.dataset{setInx}.p = {};
                    handles.dataset{setInx}.q = {};
                end

                setInx = setInx +1;
            end
        end
        
        clear('data','filenames','sequence');
    end
    set(handles.text_message,'String','Dataset creation now completed.')
end

%Copy the dataset lists to the Single Shot popup list
set(handles.popup_curDataset,'String',handles.datatitles);
set(handles.popup_curDataset,'Value',1);

%Perhaps should check to make sure all of the sequences are the same
if length(unique(seqCheck)) > 1
    warndlg('Not all datasets appear to have the same syllable sequence.');
    uiwait;
end

guidata(hObject, handles);

function listbox_datasets_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox_datasets_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single Shot Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Dataset Selection 
function popup_curDataset_Callback(hObject, eventdata, handles)
%This function kicks off the single shot analysis process. Upon selection,
%the dataset analyzed for all required values and outputs are plotted to
%the screen. 

%Clear the current axes
handles = clearSingleShotAxes(handles);

%Get the selected dataset index
ind = get(handles.popup_curDataset,'Value');

%Process the selected dataset (all computations ahould be done here). To
%save processing time, check to see if this dataset has already been
%analyzed. If it has, skip this step and go right to the plotting.
if isfield (handles, 'procCount')
    if ~(handles.procCount==ind)
        [handles.procSSData{ind}, handles.pitchSnipParams] = processDataset(handles,ind);
        handles.procCount = sort([handles.procCount; ind]);
    end
else
    [handles.procSSData{ind}, handles.pitchSnipParams] = processDataset(handles,ind);
    handles.procCount = ind;
end

%Populate the rendition spectrogram elements; Set GUI controls
rendNum = length(handles.dataset{ind}.audio);
set(handles.slider_rendCntr,'Min',1,'Max',rendNum,'Value',1,'SliderStep',[(1/rendNum), (1/rendNum)]);
set(handles.edit_endNum,'String','1')

%Plot the first rendition spectrogram
axes(handles.axes_rendSpect)
displaySpecgramQuick(handles.dataset{ind}.audio{1}, handles.fs,[0,10000],[],0);
set(handles.axes_rendSpect,'XTick',[],'YTick',[]);
xlabel(''); ylabel('');

%Plot the Stack Plot
axes(handles.axes_stackPlot)
imagesc(handles.dataset{ind}.audioCube)
set(handles.axes_stackPlot,'XTick',[],'YTick',[]);
xlabel(''); ylabel('');

%Plot Time-of-Day Histogram
axes(handles.axes_timeDistr)
h = updateTimeDistr(handles,handles.procSSData{ind}.hourVect,handles.procSSData{ind}.minVect);
xlim([0 24]);
set(h,'facecolor',[.5 .5 .5])
box off;
set(handles.axes_timeDistr,'XTick',0:4:24,'TickDir','out')
xlabel('Time of day');
ylabel('File Count')

%Populate the Alignment Paths plots
axes(handles.axes_alignPaths)
h = updatePathsPlot(handles.dataset{ind}.p,handles.dataset{ind}.q);
box off; axis tight
ylim([-60, 60])
set(handles.axes_alignPaths,'XTick',0:200:max(handles.dataset{ind}.p{1}),'YTick',-30:30:30,'TickDir','out');
xlabel('Template Time');
ylabel('<-Comp - Str->')

%Plot selected interval duration distributions
axes(handles.axes_SylGapDistr)
cla; hold on
handles = SylGapDistr(handles,ind);
hold off; box off

%Figure out what the user has asked to display
featNames =  fieldnames(handles.procSSData{ind}.syl_features{1});
k = strfind(featNames,'mean_');
parfor i = 1:length(k)
    bools(i) = ~isempty(k{i});
end
meanFeats = featNames(bools);
feat = meanFeats(get(handles.popup_featStats,'Value'));

%Update featStat plot given the user selections
axes(handles.axes_featStats)
cla; hold on
[featMat,featMeans,featStd] = updatefeatStats(handles.procSSData{ind}.syl_features,feat,handles.dataset{ind}.sLabels);
hold off; box off
xlim([0,length(handles.dataset{ind}.sLabels)+1])
ymin = round((min(featMat(:))-0.1*(abs(min(featMat(:)))))*100)/100; ymax = round((max(featMat(:))+0.1*(abs(max(featMat(:)))))*100)/100;
ylim([ymin,ymax]); mm = max(abs([ymin,ymax]));
set(handles.axes_featStats,'XTickLabel',handles.dataset{ind}.sLabels,'XTick',1:length(handles.dataset{ind}.sLabels),'TickDir','out')
%set(handles.axes_featStats,'YTick',round(100*(-mm:.25*mm:mm))/100,'YTickLabel',round(100*(-mm:.25*mm:mm))/100,'TickDir','out')
xlabel('Syl Type');
ylabel('XXX')

%Update the Pitch Snips plot
axes(handles.axes_pitchSnips)
cla; hold on
updatePitchSnips(handles.procSSData{ind}.pitchSnips, handles.dataset{ind}.sLabels);
hold off; box off
xlim([0,length(handles.dataset{ind}.sLabels)+1])
ymin = round((min(handles.procSSData{ind}.pitchSnips(:))-0.1*(abs(min(handles.procSSData{ind}.pitchSnips(:)))))*100)/100; ymax = round((max(handles.procSSData{ind}.pitchSnips(:))+0.1*(abs(max(handles.procSSData{ind}.pitchSnips(:)))))*100)/100;
ylim([ymin,ymax]); range = ymax-ymin;
set(handles.axes_pitchSnips, 'XTickLabel', handles.dataset{ind}.sLabels, 'XTick', 1:length(handles.dataset{ind}.sLabels))
set(handles.axes_pitchSnips, 'YTick', round(ymin:.25*range:ymax), 'YTickLabel', round(ymin:.25*range:ymax), 'TickDir', 'out')
xlabel('Syl Type');
ylabel('Pitch (Hz)')


%And copy the syllable type list to the frequency analysis popup box
set(handles.popup_freqSyl,'String', handles.dataset{ind}.sLabels)


guidata(hObject, handles);

function popup_curDataset_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


    %Rendition browser
function slider_rendCntr_Callback(hObject, eventdata, handles)
%Get the selected index
ind = get(handles.popup_curDataset,'Value');

%Populate the rendition spectrogram elements
axes(handles.axes_rendSpect)
rendNum = max(floor(get(handles.slider_rendCntr,'Value')),1);
set(handles.edit_endNum,'String',num2str(rendNum))
displaySpecgramQuick(handles.dataset{ind}.audio{rendNum}, handles.fs,[0,10000],[],0);
set(handles.axes_rendSpect,'XTick',[],'YTick',[]);
xlabel(''); ylabel('');

guidata(hObject, handles);

function slider_rendCntr_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function edit_endNum_Callback(hObject, eventdata, handles)
%Get the selected index
ind = get(handles.popup_curDataset,'Value');
rendMax = length(handles.dataset{ind}.audio);

%Populate the rendition spectrogram elements
axes(handles.axes_rendSpect)
rendNum = max(floor(str2num(get(handles.edit_endNum,'String'))),1);
rendNum = min(rendNum,rendMax);
set(handles.slider_rendCntr,'Value',rendNum)
set(handles.edit_endNum,'String',num2str(rendNum))
displaySpecgramQuick(handles.dataset{ind}.audio{rendNum}, handles.fs,[0,10000],[],0);
set(handles.axes_rendSpect,'XTick',[],'YTick',[]);
xlabel(''); ylabel('');

guidata(hObject, handles);

function edit_endNum_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


    %SAP Featire analysis
function popup_featStats_Callback(hObject, eventdata, handles)
ind = get(handles.popup_curDataset,'Value');
featNames =  fieldnames(handles.procSSData{ind}.syl_features{1});
k = strfind(featNames,'mean_');
for i = 1:length(k)
    bools(i) = ~isempty(k{i});
end
meanFeats = featNames(bools);
feat = meanFeats(get(handles.popup_featStats,'Value'));

axes(handles.axes_featStats)
cla; hold on
[featMat,featMeans,featStd] = updatefeatStats(handles.procSSData{ind}.syl_features,feat,handles.dataset{ind}.sLabels);
hold off; box off

xlim([0,length(handles.dataset{ind}.sLabels)+1])
ymin = round((min(featMat(:))-0.1*(abs(min(featMat(:)))))*100)/100; ymax = round((max(featMat(:))+0.1*(abs(max(featMat(:)))))*100)/100;
ylim([ymin,ymax]); mm = max(abs([ymin,ymax]));
set(handles.axes_featStats,'XTickLabel',handles.dataset{ind}.sLabels,'XTick',1:length(handles.dataset{ind}.sLabels))
set(handles.axes_featStats,'YTick',round(100*(-mm:.25*mm:mm))/100,'YTickLabel',round(100*(-mm:.25*mm:mm))/100)
xlabel('Int Type');
ylabel('XXX')

guidata(hObject, handles);

function popup_featStats_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function check_calcSAP_Callback(hObject, eventdata, handles)


    %Syl/Gap Duration analysis
function check_durSylDistr_Callback(hObject, eventdata, handles)
%Checking this box, sets the paired box to false
set(handles.check_durBySyllGap,'Value',0);
set(handles.check_durBySyllNGap,'Value',0);

%Update plot selected distributions
axes(handles.axes_SylGapDistr)
cla; hold on
[handles] = SylGapDistr(handles,get(handles.popup_curDataset,'Value'));
hold off; box off

guidata(hObject, handles);

function check_durGapDistr_Callback(hObject, eventdata, handles)
%Checking this box, sets the paired box to false
set(handles.check_durBySyllGap,'Value',0);
set(handles.check_durBySyllNGap,'Value',0);

%Update plot selected distributions
axes(handles.axes_SylGapDistr)
cla; hold on
[handles] = SylGapDistr(handles,get(handles.popup_curDataset,'Value'));
hold off; box off

guidata(hObject, handles);

function check_durAllDistr_Callback(hObject, eventdata, handles)
%Checking this box, sets the paired box to false
set(handles.check_durBySyllGap,'Value',0);
set(handles.check_durBySyllNGap,'Value',0);

%Update plot selected distributions
axes(handles.axes_SylGapDistr)
cla; hold on
[handles] = SylGapDistr(handles,get(handles.popup_curDataset,'Value'));
hold off; box off

guidata(hObject, handles);

function check_durBySyllGap_Callback(hObject, eventdata, handles)
%Checking this box, sets the paired box to false
set(handles.check_durSylDistr,'Value',0);
set(handles.check_durGapDistr,'Value',0);
set(handles.check_durAllDistr,'Value',0);
set(handles.check_durBySyllNGap,'Value',0);

%You may not uncheck the box
if ~get(handles.check_durBySyllGap,'Value')
    set(handles.check_durBySyllGap,'Value',1);
end

%Update plot selected distributions
axes(handles.axes_SylGapDistr)
cla; hold on
[handles] = SylGapDistr(handles,get(handles.popup_curDataset,'Value'));
hold off; box off

guidata(hObject, handles);

function check_durBySyllNGap_Callback(hObject, eventdata, handles)
%Checking this box, sets the paired box to false
set(handles.check_durSylDistr,'Value',0);
set(handles.check_durGapDistr,'Value',0);
set(handles.check_durAllDistr,'Value',0);
set(handles.check_durBySyllGap,'Value',0);

%You may not uncheck the box
if ~get(handles.check_durBySyllNGap,'Value')
    set(handles.check_durBySyllNGap,'Value',1);
end

%Update plot selected distributions
axes(handles.axes_SylGapDistr)
cla; hold on
[handles] = SylGapDistr(handles,get(handles.popup_curDataset,'Value'));
hold off; box off

guidata(hObject, handles);


    %All Syls freq analysis
function push_pitchParam_Callback(hObject, eventdata, handles)
%Function pops out a dialog window to gather parameters for pitch snips function
ind = get(handles.popup_curDataset,'Value');
if isfield(handles,'pitchSnipParams')
    handles.pitchSnipParams = getPitchSnipInfo(handles.pitchSnipParams, handles.dataset{ind}.sLabels);
else
    handles.pitchSnipParams = getPitchSnipInfo([], handles.dataset{ind}.sLabels);
end

guidata(hObject, handles);

function push_pitchUpdate_Callback(hObject, eventdata, handles)
%This function updates the pitchSnips plot with the currently selected values

%Collect the dataset index number and re-process the pitch values
ind = get(handles.popup_curDataset,'Value');
handles.procSSData{ind}.pitchSnips = pitchSnip(handles.dataset{ind}.audio, handles.dataset{ind}.p, handles.dataset{ind}.q, handles.dataset{ind}.templatesyllBreaks, handles.pitchSnipParams);

%Update the Pitch Snips plot
axes(handles.axes_pitchSnips)
cla; hold on
updatePitchSnips(handles.procSSData{ind}.pitchSnips, handles.dataset{ind}.sLabels);
hold off; box off
xlim([0,length(handles.dataset{ind}.sLabels)+1])
ymin = round((min(handles.procSSData{ind}.pitchSnips(:))-0.1*(abs(min(handles.procSSData{ind}.pitchSnips(:)))))*100)/100; ymax = round((max(handles.procSSData{ind}.pitchSnips(:))+0.1*(abs(max(handles.procSSData{ind}.pitchSnips(:)))))*100)/100;
ylim([ymin,ymax]); range = ymax-ymin;
set(handles.axes_pitchSnips, 'XTickLabel', handles.dataset{ind}.sLabels, 'XTick', 1:length(handles.dataset{ind}.sLabels))
set(handles.axes_pitchSnips, 'YTick', round(ymin:.25*range:ymax), 'YTickLabel', round(ymin:.25*range:ymax), 'TickDir', 'out')
xlabel('Syl Type');
ylabel('Pitch (Hz)')

guidata(hObject, handles);

    %Single Syl frequency analysis
function popup_freqSyl_Callback(hObject, eventdata, handles)
%Check to see if the current dataset has been processed; if not throw a message
ind = get(handles.popup_curDataset,'Value');
if ~isfield(handles, 'procSSData') || isempty(handles.procSSData{ind})
    set(handles.text_message, 'String', 'You must complete the single shot analysis for this dataset before selecting a syllable.')
    
else %If everything has properly loaded, continue on...
    %Get the selected syllable type
    sylTypeInd = get(handles.popup_freqSyl,'Value');
    centerFreq = str2double(get(handles.edit_centerFreq, 'String'));
    buffer = 10; %in ms

    %Calulate the pitch traces for the targetted syllable
    pitchCell = pitchTraces(handles.dataset{ind}.audio, handles.dataset{ind}.p, handles.dataset{ind}.q, handles.dataset{ind}.templatesyllBreaks, sylTypeInd, centerFreq, buffer);
    
    %Convert the cell to matrix, padding the end with NaNs as needed
    maxLength = max(cellfun(@(x)numel(x),pitchCell));
    pitchMatrix = cell2mat(cellfun(@(x)cat(2,x,NaN(1,maxLength-length(x))),pitchCell,'UniformOutput',false));
    
    %Calc stats on the pitch traces
    pitchMean = nanmean(pitchMatrix,1);
    pitchCV = nanstd(pitchMatrix,1)./pitchMean;
    
    %Plot the output on to the axes
    axes(handles.axes_Freq)
    cla; hold on; box off;
    xs = (-buffer+1):1:(maxLength-buffer);
    plot(xs', pitchMatrix,'b','LineWidth', 0.25)
    plot(xs, pitchMean, 'k', 'LineWidth', 2.0)
    
    set(handles.axes_Freq, 'TickDir', 'out')
    xlabel('Time from Onset (ms)');
    ylabel('Pitch (Hz)')
    hold off
    
end

guidata(hObject, handles);




function popup_freqSyl_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end    
    
function slider_centerFreq_Callback(hObject, eventdata, handles)

function slider_centerFreq_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function edit_centerFreq_Callback(hObject, eventdata, handles)

function edit_centerFreq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function slider_winOffset_Callback(hObject, eventdata, handles)

function slider_winOffset_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function edit_winOffset_Callback(hObject, eventdata, handles)

function edit_winOffset_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Longitudinal Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function push_longProcess_Callback(hObject, eventdata, handles)
%Clear all longitudinal axes
[handles] = clearLongitudinalAxes(handles);

%Calculate reference datenum
handles.long.refDateNum = datenum(get(handles.edit_startDate,'String'),'yymmdd');

%Cycle through each dataset and calulate the required quantities; skip
%those that have already been analyzed. Locate number of datasets
% if iscell(handles.datasetFilename)
%     numDatasets = size(handles.datasetFilename,2);
% else
%     numDatasets = size(handles.datasetFilename,1);
% end

numDatasets = size(handles.datatitles,1);

for ind=1:numDatasets
    %All datasets in the longitudinal analysis need to be based on the same sequence, so compare
    %each sequence to the first dataset and throw an error if they don't match
    if ~isfield(handles.long, 'sequence')
        handles.long.sequence = handles.dataset{ind}.sequence;
    else
        if handles.long.sequence ~= handles.dataset{ind}.sequence
            set(handles.text_message,'String',['The sequence for dataset ' handles.datatitles{ind} ...
                ' does not match previous sequences. Check the file and try again.']);
            return
        end
    end
    
    %Process the selected dataset. If it has, skip this step.
    if isfield (handles, 'procCount')
        if ~(handles.procCount==ind)
            handles.procSSData{ind} = processDataset(handles,ind);
            handles.procCount = sort([handles.procCount; ind]);
        end
    else
        handles.procSSData{ind} = processDataset(handles,ind);
        handles.procCount = ind;
    end

    %Temporary for trouble shooting.
    %handles.procSSData{ind} = processDataset(handles,ind);

    %Calculate the number of days from the reference
    handles.procSSData{ind}.recDateNum = handles.dataset{ind}.datenum - handles.long.refDateNum;
end

%Capture the syl-gap selections and copy labels to the popups.
if get(handles.check_longSylGap,'Value')
    contents = ['NA'; handles.dataset{1}.sgLabels; 'Total'];
    set(handles.popup_temporalBlue,'String',contents)
    set(handles.popup_temporalRed,'String',contents)
    set(handles.popup_temporalGreen,'String',contents)
elseif get(handles.check_longSylNGap,'Value')
    contents = ['NA'; handles.dataset{1}.sNgLabels; 'Total'];
    set(handles.popup_temporalBlue,'String',contents)
    set(handles.popup_temporalRed,'String',contents)
    set(handles.popup_temporalGreen,'String',contents)
end

%Copy labels to spectral popups
contents = ['NA'; handles.dataset{1}.sLabels];
set(handles.popup_spectralBlue,'String',contents)
set(handles.popup_spectralRed,'String',contents)
set(handles.popup_spectralGreen,'String',contents)

contents = ['Pitch Snip'; 'SAP Var'; handles.featureSet];
set(handles.popup_featSelect,'String',contents)

%Set output message when processing is finished
set(handles.text_message,'String','Completed processing all datasets. Select options to display output.')

guidata(hObject, handles);

function edit_startDate_Callback(hObject, eventdata, handles)

function edit_startDate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
                    %Temporal
function check_longSylGap_Callback(hObject, eventdata, handles)
%Checking this box, sets the paired box to false
set(handles.check_longSylNGap,'Value',0);

%You may not uncheck the box
if ~get(handles.check_longSylGap,'Value')
    set(handles.check_longSylGap,'Value',1);
end

%Update plot selected distributions
contents = ['NA'; handles.dataset{1}.sgLabels; 'Total'];
set(handles.popup_temporalBlue,'String',contents,'Value',1)
set(handles.popup_temporalRed,'String',contents,'Value',1)
set(handles.popup_temporalGreen,'String',contents,'Value',1)
    
guidata(hObject, handles);

function check_longSylNGap_Callback(hObject, eventdata, handles)
%Checking this box, sets the paired box to false
set(handles.check_longSylGap,'Value',0);

%You may not uncheck the box
if ~get(handles.check_longSylNGap,'Value')
    set(handles.check_longSylNGap,'Value',1);
end

%Update plot selected distributions
contents = ['NA'; handles.dataset{1}.sNgLabels; 'Total'];
set(handles.popup_temporalBlue,'String',contents,'Value',1)
set(handles.popup_temporalRed,'String',contents,'Value',1)
set(handles.popup_temporalGreen,'String',contents,'Value',1)

guidata(hObject, handles);

function popup_temporalBlue_Callback(hObject, eventdata, handles)
%Call function to update the longitudinal plot
handles = updateLongTemporal(handles);

guidata(hObject, handles);

function popup_temporalBlue_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_temporalRed_Callback(hObject, eventdata, handles)
%Call function to update the longitudinal plot
handles = updateLongTemporal(handles);

guidata(hObject, handles);

function popup_temporalRed_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_temporalGreen_Callback(hObject, eventdata, handles)
%Call function to update the longitudinal plot
handles = updateLongTemporal(handles);

guidata(hObject, handles);

function popup_temporalGreen_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function check_tempNormal_Callback(hObject, eventdata, handles)
%Call function to update the longitudinal plot
handles = updateLongTemporal(handles);

guidata(hObject, handles);

function check_tempCV_Callback(hObject, eventdata, handles)
%Call function to update the longitudinal plot
handles = updateLongTemporal(handles);

guidata(hObject, handles);

                    %Spectral
function popup_featSelect_Callback(hObject, eventdata, handles)
%Call function to update the longitudinal plot
handles = updateLongSpectral(handles);

guidata(hObject, handles);

function popup_featSelect_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_spectralBlue_Callback(hObject, eventdata, handles)
%Call function to update the longitudinal plot
handles = updateLongSpectral(handles);

guidata(hObject, handles);

function popup_spectralBlue_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_spectralRed_Callback(hObject, eventdata, handles)
%Call function to update the longitudinal plot
handles = updateLongSpectral(handles);

guidata(hObject, handles);

function popup_spectralRed_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_spectralGreen_Callback(hObject, eventdata, handles)
%Call function to update the longitudinal plot
handles = updateLongSpectral(handles);

guidata(hObject, handles);

function popup_spectralGreen_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

                    %Optional
function popup_optionalBlue_Callback(hObject, eventdata, handles)

function popup_optionalBlue_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_optionalRed_Callback(hObject, eventdata, handles)

function popup_optionalRed_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_optionalGreen_Callback(hObject, eventdata, handles)

function popup_optionalGreen_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function check_optionalNormal_Callback(hObject, eventdata, handles)


function check_specNormal_Callback(hObject, eventdata, handles)
%Call function to update the longitudinal plot
handles = updateLongSpectral(handles);

guidata(hObject, handles);

function check_specCV_Callback(hObject, eventdata, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_exportRaw_Callback(hObject, eventdata, handles)
%Generate a new figure and export the data from the longitudinal plots
% figure;
% subplot(2,1,1)
% updateLongTemporal(handles);
% 
% subplot(2,1,2)
% updateLongSpectral(handles);

function push_exportProc_Callback(hObject, eventdata, handles)

data.procSSData = handles.procSSData;

if iscell(handles.datasetFilename)
    numDatasets = size(handles.datasetFilename,2);
else
    numDatasets = size(handles.datasetFilename,1);
end

for ind=1:numDatasets
    data.day(ind) = handles.procSSData{ind}.recDateNum;
end
data.filenames = handles.datasetFilename;
data.pitchParams = handles.pitchSnipParams;

%Requst location and save to file.
[SaveName,SavePath] = uiputfile('*.mat');
save([SavePath SaveName],'data');

guidata(hObject, handles);

function push_exportLong_Callback(hObject, eventdata, handles)
%Parse the data for the current display
%Cycle through each dataset and calulate the required quantities; skip
%those that have already been analyzed. Locate number of datasets
data = handles.long;

if iscell(handles.datasetFilename)
    numDatasets = size(handles.datasetFilename,2);
else
    numDatasets = size(handles.datasetFilename,1);
end

for ind=1:numDatasets
    data.day(ind) = handles.procSSData{ind}.recDateNum;
end
data.filenames = handles.datasetFilename;
data.pitchParams = handles.pitchSnipParams;

%Requst location and save to file.
[SaveName,SavePath] = uiputfile('*.mat');
save([SavePath SaveName],'data');

guidata(hObject, handles);

function push_screenCapture_Callback(hObject, eventdata, handles)



function check_loadFull_Callback(hObject, eventdata, handles)
%Turn these gui objects on/off as they are/aren't needed
if get(handles.check_loadFull, 'Value')
    set(handles.edit_parseTimes, 'Visible', 'off')
    set(handles.edit_numMotifs, 'Visible', 'off')
    set(handles.text23, 'Visible', 'off')
    set(handles.text24, 'Visible', 'off')
else
    set(handles.edit_parseTimes, 'Visible', 'on')
    set(handles.edit_numMotifs, 'Visible', 'on')
    set(handles.text23, 'Visible', 'on')
    set(handles.text24, 'Visible', 'on')
end

guidata(hObject, handles);


function edit_parseTimes_Callback(hObject, eventdata, handles)

function edit_parseTimes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_numMotifs_Callback(hObject, eventdata, handles)

function edit_numMotifs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
