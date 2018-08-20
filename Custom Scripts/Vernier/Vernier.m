function varargout = Vernier(varargin)
% VERNIER MATLAB code for Vernier.fig
%      VERNIER, by itself, creates a new VERNIER or raises the existing
%      singleton*.
%
%      H = VERNIER returns the handle to a new VERNIER or the handle to
%      the existing singleton*.
%
%      VERNIER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VERNIER.M with the given input arguments.
%
%      VERNIER('Property','Value',...) creates a new VERNIER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Vernier_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Vernier_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Vernier

% Last Modified by GUIDE v2.5 27-Mar-2013 12:58:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Vernier_OpeningFcn, ...
                   'gui_OutputFcn',  @Vernier_OutputFcn, ...
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


% --- Executes just before Vernier is made visible.
function Vernier_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Vernier (see VARARGIN)

% Choose default command line output for Vernier
handles.output = hObject;

%Clear the crap off the axes
handles = clearAxes(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Vernier wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Vernier_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stand Alone Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = clearAxes(handles)
%Clear all of the GUI axes
axes(handles.axes_breaks)
cla
set(gca,'XTick',[],'YTick',[],'TickDir','out','Box','off')
xlabel(''); ylabel('')

axes(handles.axes_spec1)
cla
set(gca,'XTick',[],'YTick',[],'TickDir','out','Box','off')
xlabel(''); ylabel('')

axes(handles.axes_spec2)
cla
set(gca,'XTick',[],'YTick',[],'TickDir','out','Box','off')
xlabel(''); ylabel('')

axes(handles.axes_audio1)
cla
set(gca,'XTick',[],'YTick',[],'TickDir','out','Box','off')
xlabel(''); ylabel('')

axes(handles.axes_audio2)
cla
set(gca,'XTick',[],'YTick',[],'TickDir','out','Box','off')
xlabel(''); ylabel('')

axes(handles.axes_audio3)
cla
set(gca,'XTick',[],'YTick',[],'TickDir','out','Box','off')
xlabel(''); ylabel('')

axes(handles.axes_neuro1)
cla
set(gca,'XTick',[],'YTick',[],'TickDir','out','Box','off')
xlabel(''); ylabel('')

axes(handles.axes_neuro2)
cla
set(gca,'XTick',[],'YTick',[],'TickDir','out','Box','off')
xlabel(''); ylabel('')

axes(handles.axes_neuro3)
cla
set(gca,'XTick',[],'YTick',[],'TickDir','out','Box','off')
xlabel(''); ylabel('')

function data = createSummary(data)
%Preallocate required variables
intervals = [];
alignedAudio = [];
HRAudio = [];
alignedNeuro = [];

%Work through each rendition to extract summary data
specWin = 220; %5ms
specAdv = 44; %1ms
fs = 44150;
for i = 1:length(data.audio)
    %Retrieve the rendition warping path
    path = round([data.p{i},data.q{i}]);
    
    %Calculate all interval lengths for the rendition
    syllBreaks = getWarpedStarts(path,data.templateSyllBreaks);
    intervals = [intervals;diff(syllBreaks)'];
    
    %Generate rendition spectrogram
    [~,~,~,P] = spectrogram((data.audio{i}/(sqrt(mean(data.audio{i}.^2)))),specWin,specWin-specAdv,512,fs);
    rendLength = size(P,2);
    
    %Generate high resolution spectrogram
    [~,~,~,PHR] = spectrogram((data.audio{i}/(sqrt(mean(data.audio{i}.^2)))),specWin,specWin-specAdv,4096,fs);
    
    %Create sparse path from interval lengths 
    path_t = [1; data.templateSyllBreaks; data.templateLength];
    path_r = [1; syllBreaks; rendLength];
    sp_path = [path_t,path_r];
    
    %Align
    alignedAudio(i,:,:) = alignSeriesSTW(abs(P(4:94,:)),sp_path);
    
    HRAudio(i,:,:) = alignSeriesSTW(abs(PHR(1:1000,:)),sp_path);
    
    if isfield(data,'neuro')
        %Calc Neural power and smooth with boxcar
        pow = mean(windowTS_spec(data.neuro{i}.^2,specWin,specAdv));
      
        %Warp the rendition neural power using the sparse path
        alignedNeuro(i,:) = alignSeriesSTW(pow,sp_path);
    end
    
end
%Log-transform and normalize mean spectrogram
t = 10*log10(squeeze(mean(alignedAudio,1)));
meanSpec = t./mean(t(:));

t = 10*log10(squeeze(mean(HRAudio,1)));
meanHRSpec = t./mean(t(:));

meanIntervals = mean(intervals,1);

%Get similarly aligned/warped neural trace
if isfield(data,'neuro')
    meanNeuro=squeeze(mean(alignedNeuro,1));
end

%Copy useful stuff to data structure
data.meanSpec = meanSpec;
data.meanHRSpec = meanHRSpec;
data.meanIntervals = meanIntervals;
if isfield(data,'neuro')
    data.meanNeuro=meanNeuro;
end

function data = unwarpSeries(data)
%Create the inverse path
meanPath = round([data.templateSyllBreaks(1), cumsum(data.meanIntervals)+data.templateSyllBreaks(1)])';
invPath = [data.templateSyllBreaks, meanPath];
data.invPath = fliplr([1,1;invPath;data.templateLength,size(data.meanSpec,2)]);

%Unwarp the series/spectrogram with the sparse inverse path
data.un_meanSpec = alignSeriesSTW(data.meanSpec,data.invPath);
data.un_meanHRSpec = alignSeriesSTW(data.meanHRSpec,data.invPath);
if isfield(data,'neuro')
    data.un_meanNeuro = alignSeriesSTW(data.meanNeuro,data.invPath);
end

function edges = getNewEdges(Spec,thresh,starts)
%Collapse into power timeseries
specPowTS = -10*log10(sum(Spec));

%Use EM to fit 2 gaussians to distribution
gMix = gmdistribution.fit(specPowTS',2);
[powerM,pnt] = sort(gMix.mu);
temp = sqrt(squeeze(gMix.Sigma));
powerStd = temp(pnt);
thresh = powerM(1)+thresh*powerStd(1); %thresh is SDs above the silent Gaussian mean

%Find the crossing points to estimate the onset and offsets of
%motifs
ind = 2:length(specPowTS);
Crossings.Up = find(specPowTS(ind)>thresh & specPowTS(ind-1)<thresh);
Crossings.Down = find(specPowTS(ind)<thresh & specPowTS(ind-1)>thresh);

onStarts = starts(1:2:end);
offStarts = starts(2:2:end);

edges = [];
 for k = 1:length(onStarts)
     [offset, pntr] = min(abs(Crossings.Up-onStarts(k)));
     edges = [edges;Crossings.Up(pntr)];
     
     [offset, pntr] = min(abs(Crossings.Down-offStarts(k)));
     edges = [edges;Crossings.Down(pntr)];
 end

function [powRatio1,powRatio2] = harmonicPower(handles)
%Retrive parameters for calcs
pitchCenter = str2double(get(handles.edit_pitchCenter,'String')); %mean pitch for target syllable
numHarm = str2double(get(handles.edit_numHarm,'String')); %number of harmonics to include
specWin = 220; %5ms
specAdv = 44; %1ms
fs = 44150;

%Get the frequency bins for the high-res spectrograms
[~,freqBins,~,~] = spectrogram(handles.data1.audio{1},specWin,specWin-specAdv,4096,fs);

%Set up indices for the target pitch harmonics
indx = [];
for  i = 1:numHarm
    diffVect = abs(freqBins-(i*pitchCenter));
    [~,pntr] = min(diffVect);
    indx = [indx;pntr];
end
offset = round(mean(diff(indx)/2));
indx_off = indx+offset;

%Harmonic power analysis
PowerON1 =  handles.data1.un_meanHRSpec(indx,:); 
PowerOFF1 =  handles.data1.un_meanHRSpec(indx_off,:); 
PowerON2 =  handles.data2.un_meanHRSpec(indx,:); 
PowerOFF2 =  handles.data2.un_meanHRSpec(indx_off,:); 

temp1 = [];
temp2 = [];
for i=1:numHarm
    temp1(i,:) = PowerON1(i,:)./PowerOFF1(i,:);
    temp2(i,:) = PowerON2(i,:)./PowerOFF2(i,:);
end
smWin = str2double(get(handles.edit_smWindow,'String'));
if smWin>0
    powRatio1 = smooth(mean(temp1,1),smWin);
    powRatio2 = smooth(mean(temp2,1),smWin);
else
    powRatio1 = mean(temp1,1);
    powRatio2 = mean(temp2,1);
end

function [truncIndx1,truncIndx2] = truncSeries(anchors)
posLength1 = anchors(3,1)-anchors(2,1);
posLength2 = anchors(3,2)-anchors(2,2);
posTrunc = min(posLength1,posLength2);

negLength1 = anchors(2,1)-anchors(1,1);
negLength2 = anchors(2,2)-anchors(1,2);
negTrunc = min(negLength1,negLength2);

truncIndx1 = anchors(2,1)-negTrunc:anchors(2,1)+posTrunc;
truncIndx2 = anchors(2,2)-negTrunc:anchors(2,2)+posTrunc;

function pCorr = winPearson(ts1,ts2,window)
%Window each timeseries
[ts1_mat] = windowTS(ts1,window,1,'pad','boxcar');
[ts2_mat] = windowTS(ts2,window,1,'pad','boxcar');

%Work through each matching window and calc Pearon correlation
winNum = size(ts1_mat,1);
for i=1:winNum-1
    ind = ~isnan(ts1_mat(i,:));
    pCorr(i) = (corr(ts1_mat(i,ind)',ts2_mat(i,ind)','type','Pearson'));
end

function [pnts,locs] = salientPoints(TS,minSpaceDist,minGradient,checkDist)
%Takes in a timeseries (TS) and return a set of points (in x-y space) that gives
%the salient peaks and valleys.  This method assumes that all values of TS
%are positive.

lengthTS = length(TS);
% minSpaceDist = 3;
% minGradient = 0.005;
% checkDist = 3;

%Get Peaks
[peaks,pLocs]=findpeaks(TS,'MINPEAKDISTANCE',minSpaceDist);
pmask = zeros(length(peaks),1);
for i = 1:length(peaks)
    steepLeft = (TS(max(1,pLocs(i)-checkDist)) <= (peaks(i)-minGradient));
    steepRight = (TS(min(pLocs(i)+checkDist,lengthTS))) <= (peaks(i)-minGradient);
    if steepLeft && steepRight
        pmask(i) = 1;
    end
end

%Get Valleys
[valleys,vLocs]=findpeaks(-TS,'MINPEAKDISTANCE',minSpaceDist);
valleys = -valleys;
vmask = zeros(length(valleys),1);
for i = 1:length(valleys)
    steepLeft = (TS(max(1,vLocs(i)-checkDist)) >= (valleys(i)+minGradient));
    steepRight = (TS(min(vLocs(i)+checkDist,lengthTS))) >= (valleys(i)+minGradient);
    if steepLeft && steepRight
        vmask(i) = 1;
    end
end

%Merge peaks and valleys into single set
tempLoc = [pLocs(logical(pmask)),vLocs(logical(vmask))];
tempPnts = [peaks(logical(pmask)),valleys(logical(vmask))];

[locs,indx] = sort(tempLoc);
pnts = tempPnts(indx);

%Plot output...
figure
plot(TS); hold on; axis tight; ylim([0,1])
scatter(pLocs,peaks,'or')
scatter(vLocs,valleys,'og')
scatter(locs,pnts,'+k')

function handles = plotAlignments(handles)
%Calculate the shift indices required to align plots tot he selected intervals
handles.data1.timeVect = 1:handles.data.d2dAudioPath(end,1);
handles.data2.timeVect = 1:handles.data.d2dAudioPath(end,2);
targetTimes = [handles.data1.sylEdges(handles.data.targetInt(1)),handles.data2.sylEdges(handles.data.targetInt(1))];

%Set plotting limits
negLim = -1*max(targetTimes);
posLim = max(handles.data.d2dAudioPath(end,1)-targetTimes(1),handles.data.d2dAudioPath(end,2)-targetTimes(2));

%Show interval break points on the top axis
axes(handles.axes_breaks)
cla
hold on
for i = 1:length(handles.data1.sylEdges)
    text(handles.data1.sylEdges(i)-targetTimes(1)-8, mean(ylim), num2str(i), 'FontWeight', 'bold', 'FontSize', 16);
end
xlim([negLim,posLim])
set(gca,'XTick',[],'YTick',[],'TickDir','out');
set(gca,'Box','off')
hold off

%Plot aligned Spectrograms
axes(handles.axes_spec1)
cla; hold on
imagesc(-targetTimes(1),-0.5,-handles.data1.un_meanSpec); axis xy
barlims = ones(length(handles.data1.sylEdges),1)*ylim;
line([handles.data1.sylEdges-targetTimes(1),handles.data1.sylEdges-targetTimes(1)]',barlims','Color','k','LineStyle','--','LineWidth',1.5)
axis tight
xlim([negLim,posLim])
set(gca,'XTick',[],'YTick',[],'TickDir','out','Box','off');

axes(handles.axes_spec2)
cla; hold on
imagesc(-targetTimes(2),-0.5,-handles.data2.un_meanSpec); axis xy
line([handles.data2.sylEdges-targetTimes(2),handles.data2.sylEdges-targetTimes(2)]',barlims','Color','k','LineStyle','--','LineWidth',1.5)
axis tight
xlim([negLim,posLim])
xticks = unique([fliplr(-1*(0:100:-negLim)),0:100:posLim]);
set(gca,'XTick',xticks,'YTick',[],'TickDir','out','Box','off');

%Plot aligned audio power envelops
axes(handles.axes_audio1)
cla; hold on
plot(handles.data1.timeVect-targetTimes(1),mat2gray(-10*log10(sum(handles.data1.un_meanSpec,1))))
plot(handles.data2.timeVect-targetTimes(2),mat2gray(-10*log10(sum(handles.data2.un_meanSpec,1))),'r')
plot(handles.data1.timeVect-targetTimes(1),mat2gray(-10*log10(sum(handles.data2.re_meanSpec,1))),'g')
line([0,0],[0,1],'Color',[0 0 0],'LineStyle','--','LineWidth',1.5)
axis tight;
xlim([negLim,posLim])
set(gca,'XTick',[],'YTick',[0.5,1],'TickDir','out','Box','off')

%Plot aligned audio power correlations
axes(handles.axes_audio2)
cla; hold on
plot(handles.data1.timeVect(handles.data1.trncInd)-targetTimes(1),handles.data.un_Acorr,'k')
plot(handles.data1.timeVect(handles.data1.trncInd)-targetTimes(1),handles.data.re_Acorr,'g')
line([0,0],[-1,1],'Color',[0 0 0],'LineStyle','--','LineWidth',1.5)
axis tight;
xlim([negLim,posLim])
ylim([-1,1.1])
set(gca,'XTick',xticks,'YTick',[-1,0,1],'TickDir','out','Box','off')
hold off

%Plot aligned pitch/harmonic power ratios
axes(handles.axes_audio3)
cla; hold on
if get(handles.check_harmPow,'Value')
    plot(handles.data1.timeVect-targetTimes(1),handles.data1.harmPowRatio)
    plot(handles.data2.timeVect-targetTimes(2),handles.data2.harmPowRatio,'r')
    line([0,0],[0.8,1.2],'Color',[0 0 0],'LineStyle','--','LineWidth',1.5)
    xlim([negLim,posLim])
    ylim([.8,1.2])
    set(gca,'XTick',[],'YTick',[0.85,1,1.15],'TickDir','out','Box','off')
elseif get(handles.check_harmDiff,'Value')
    diffplot = (handles.data1.harmPowRatio(handles.data1.trncInd)-handles.data2.harmPowRatio(handles.data2.trncInd)).^2;
    area(handles.data1.timeVect(handles.data1.trncInd)-targetTimes(1),diffplot);
    axis tight
    line([0,0],ylim,'Color',[0 0 0],'LineStyle','--','LineWidth',1.5)
    xlim([negLim,posLim])
    set(gca,'XTick',[],'YTick',ylim,'TickDir','out','Box','off')
elseif get(handles.check_pitch,'Value')
    
end

%Plot aligned neural power envelops
axes(handles.axes_neuro1)
cla;
if isfield(handles.data2,'neuro')
    hold on
    maxmax = max(max(handles.data1.un_meanNeuro./mean(handles.data1.un_meanNeuro)),max(handles.data2.un_meanNeuro./mean(handles.data2.un_meanNeuro)));
    minmin = min(min(handles.data1.un_meanNeuro./mean(handles.data1.un_meanNeuro)),min(handles.data2.un_meanNeuro./mean(handles.data2.un_meanNeuro)));
    %plot(tShift1(Indx1),mat2gray(un_av_neuroTS1(Indx1)./mean(un_av_neuroTS1(Indx1)),[minmin,maxmax]))
    plot(handles.data1.timeVect-targetTimes(1),mat2gray(handles.data1.un_meanNeuro./mean(handles.data1.un_meanNeuro),[minmin,maxmax]))
    plot(handles.data2.timeVect-targetTimes(2),mat2gray(handles.data2.un_meanNeuro./mean(handles.data2.un_meanNeuro),[minmin,maxmax]),'r')
    plot(handles.data1.timeVect-targetTimes(1),mat2gray(handles.data2.re_meanNeuro./mean(handles.data2.re_meanNeuro),[minmin,maxmax]),'g')
    line([0,0],[0,1],'Color',[0 0 0],'LineStyle','--','LineWidth',1.5)
end
axis tight;
xlim([negLim,posLim])
set(gca,'XTick',[],'YTick',[0.5,1],'TickDir','out','Box','off')
hold off

%Plot neural power correlations
axes(handles.axes_neuro2)
cla;
if isfield(handles.data2,'neuro')
    hold on
    plot(handles.data1.timeVect(handles.data1.trncInd)-targetTimes(1),handles.data.un_Ncorr,'k')
    plot(handles.data1.timeVect(handles.data1.trncInd)-targetTimes(1),handles.data.re_Ncorr,'g')
    line([0,0],[-1,1],'Color',[0 0 0],'LineStyle','--','LineWidth',1.5)
end
axis tight;
xlim([negLim,posLim])
ylim([-1,1])
set(gca,'XTick',[],'YTick',[-1,0,1],'TickDir','out','Box','off')
hold off

%Plot 
axes(handles.axes_neuro3)
cla; hold on
if isfield(handles.data2,'neuro')
    line([handles.data1.sylEdges-targetTimes(1),handles.data1.sylEdges-targetTimes(1)]',((barlims./3)+0)','Color','b','LineStyle','--','LineWidth',1.5)
    line([handles.data2.sylEdges-targetTimes(2),handles.data2.sylEdges-targetTimes(2)]',((barlims./3)+(barlims(1,2)./3))','Color','r','LineStyle','--','LineWidth',1.5)
    line([handles.data.d2dNeuroPath(2:end-1,2)-targetTimes(1),handles.data.d2dNeuroPath(2:end-1,2)-targetTimes(1)]',((barlims./3)+(2*barlims(1,2)./3))','Color','g','LineStyle','--','LineWidth',1.5)
end
axis tight;
xlim([negLim,posLim])
set(gca,'XTick',xticks,'YTick',[0.5,1],'TickDir','out','Box','off')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_loadDay1_Callback(hObject, eventdata, handles)
%Select 'Process Dataset' (i.e., output from StretchEm)
[fname,pathLoc] = uigetfile('C:\Users\Tim\Desktop\Ali Datasets Again\*.mat');

%Load 'data' structure from the file and extract the important fields
load([pathLoc, fname],'data')
if ~exist('data','var')
    set(handles.text_message,'String',[fname ' dataset not properly formatted. Check file and retry.']);
    beep
    return
end

handles.data1.audio = data.PPaudio;
handles.data1.p = data.p;
handles.data1.q = data.q;
a = data.templatesyllBreaks'; handles.data1.templateSyllBreaks = a(:);
handles.data1.templateLength = size(data.template,2);

if isfield(data,'PPneuro')
    handles.data1.neuro = data.PPneuro;
end

%Set filename and status indicator
set(handles.text_dataset1,'String',fname);
set(handles.text_message,'String',[fname ' successfully loaded']);

%To keep memory free, delect the original 'data' structure from memore
clear data

guidata(hObject, handles);

function push_loadDay2_Callback(hObject, eventdata, handles)
%Select 'Process Dataset' (i.e., output from StretchEm)
[fname,pathLoc] = uigetfile('C:\Users\Tim\Desktop\Ali Datasets Again\*.mat');

%Load 'data' structure from the file and extract the important fields
load([pathLoc, fname],'data')
if ~exist('data','var')
    set(handles.text_message,'String',[fname ' dataset not properly formatted. Check file and retry.']);
    beep
    return
end

handles.data2.audio = data.PPaudio;
handles.data2.p = data.p;
handles.data2.q = data.q;
a = data.templatesyllBreaks'; handles.data2.templateSyllBreaks = a(:);
handles.data2.templateLength = size(data.template,2);

if isfield(data,'PPneuro')
    handles.data2.neuro = data.PPneuro;
end

%Set filename and status indicator
set(handles.text_dataset2,'String',fname);
set(handles.text_message,'String',[fname ' successfully loaded']);

%To keep memory free, delect the original 'data' structure from memore
clear data

guidata(hObject, handles);

function push_clearDay1_Callback(hObject, eventdata, handles)
%Confirm selection
button = questdlg('Are you sure you want to clear dataset?','Clear Dataset?','Clear''em!','Nooooo!','Nooooo!');

if (strcmp(button,'Clear''em!'))
    %Clear data in handles structure
    handles = rmfield(handles,'data1');

    %Set filename and status indicator
    set(handles.text_dataset1,'String','None Selected');
    set(handles.text_message,'String','Dataset cleared from memory');
else
    set(handles.text_message,'String','Dataset is unchanged.');
end

guidata(hObject, handles);

function push_clearDay2_Callback(hObject, eventdata, handles)
%Confirm selection
button = questdlg('Are you sure you want to clear dataset?','Clear Dataset?','Clear''em!','Nooooo!','Nooooo!');

if (strcmp(button,'Clear''em!'))
    %Clear data in handles structure
    handles = rmfield(handles,'data2');

    %Set filename and status indicator
    set(handles.text_dataset2,'String','None Selected');
    set(handles.text_message,'String','Dataset cleared from memory');
else
    set(handles.text_message,'String','Dataset is unchanged.');
end

guidata(hObject, handles);

function push_processData_Callback(hObject, eventdata, handles)
%Check to see that both day1 and day2 datasets are loaded
if ~isfield(handles,'data1')
    set(handles.text_message,'String','Must load Day 1 data before processing');
    beep
    return
end
if ~isfield(handles,'data2')
    set(handles.text_message,'String','Must load Day 2 data before processing');
    beep
    return
end

%Clear Axes
handles = clearAxes(handles);

%Create summary data for each day
handles.data1 = createSummary(handles.data1);
handles.data2 = createSummary(handles.data2);

%Show the unaligned spectrograms for each day
axes(handles.axes_spec1)
imagesc(-handles.data1.meanSpec)
axis xy; axis tight
set(gca,'XTick',[],'YTick',[],'TickDir','out');
set(gca,'Box','off')

axes(handles.axes_spec2)
imagesc(-handles.data2.meanSpec)
axis xy; axis tight
set(gca,'XTick',[],'YTick',[],'TickDir','out');
set(gca,'Box','off')

%Set status message
set(handles.text_message,'String','Dataset processing complete');
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alignment Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_unwarp_Callback(hObject, eventdata, handles)
%Check to see that prerequisite steps completed
if ~isfield(handles,'data1') || ~isfield(handles,'data2')
    set(handles.text_message,'String','Must load data before processing');
    beep
    return
end
if ~isfield(handles.data1,'meanSpec')
    set(handles.text_message,'String','Must load process datasets before unwarping');
    beep
    return
end

%Clear Axes
handles = clearAxes(handles);

%Unwarp all data using the mean paths
handles.data1 = unwarpSeries(handles.data1);
handles.data2 = unwarpSeries(handles.data2);

%Calculate interval edges by the chosen method
if ~get(handles.check_newEdges,'Value')
    %Simply use the syllable edges as defined by the warping paths
    handles.data1.sylEdges = handles.data1.invPath(2:end-1,1);
    handles.data2.sylEdges = handles.data2.invPath(2:end-1,1);
else
    %Recalculate syllable edges
    thresh = str2double(get(handles.edit_edgeThresh,'String'));
    handles.data1.sylEdges = getNewEdges(handles.data1.un_meanSpec,thresh,handles.data1.invPath(2:end-1,1));
    handles.data2.sylEdges = getNewEdges(handles.data2.un_meanSpec,thresh,handles.data2.invPath(2:end-1,1));
end

%Show the unaligned spectrograms for each day
axes(handles.axes_spec1)
cla
imagesc(-handles.data1.un_meanSpec)
barlims = ones(size(handles.data1.sylEdges),1)*[1,91];
line([handles.data1.sylEdges,handles.data1.sylEdges]',barlims','Color','k','LineStyle','--','LineWidth',1.5)
axis xy; axis tight
x = xlim;
set(gca,'XTick',[],'YTick',[],'TickDir','out');
set(gca,'Box','off')

axes(handles.axes_spec2)
cla
imagesc(-handles.data2.un_meanSpec)
line([handles.data2.sylEdges,handles.data2.sylEdges]',barlims','Color','k','LineStyle','--','LineWidth',1.5)
axis xy; axis tight
set(gca,'YTick',[],'TickDir','out');
set(gca,'Box','off')

%Show interval break points on the top axis
axes(handles.axes_breaks)
cla
hold on
for i = 1:length(handles.data1.sylEdges)
    text(handles.data1.sylEdges(i)-5, mean(ylim), num2str(i), 'FontWeight', 'bold', 'FontSize', 16);
end
xlim(x)
set(gca,'XTick',[],'YTick',[],'TickDir','out');
set(gca,'Box','off')
hold off

%Set status message
set(handles.text_message,'String','Data unwarped and display updated');

guidata(hObject, handles);

function push_align_Callback(hObject, eventdata, handles)
%Get target/alignment interval markers from user selections
handles.data.targetInt = [str2double(get(handles.edit_targetStart,'String')),str2double(get(handles.edit_targetEnd,'String'))];

%If selected, global stretch day 2 to match the length of non-targeted region
if get(handles.check_matchGlobal,'Value')
    %Get length of each of the non-targeted intervals
    dur1 = handles.data1.sylEdges(handles.data.targetInt(1))-handles.data1.sylEdges(1);
    dur2 = handles.data2.sylEdges(handles.data.targetInt(1))-handles.data2.sylEdges(1);
    
    %Strect factor is simply the ratio of the two
    handles.data.globalStretch = dur1/dur2;
    
    %Apply global stretch to spectrograms, neural power, and edges
    origLength = size(handles.data2.un_meanSpec,2);
    origLengthVect = 1:origLength;
    newLengthVect = linspace(1,origLength,round(origLength*handles.data.globalStretch));
    handles.data2.un_meanSpec = interp1(origLengthVect,handles.data2.un_meanSpec',newLengthVect)';
    handles.data2.un_meanHRSpec = interp1(origLengthVect,handles.data2.un_meanHRSpec',newLengthVect)';
    if isfield(handles.data2,'neuro')
        handles.data2.un_meanNeuro = interp1(origLengthVect,handles.data2.un_meanNeuro,newLengthVect);
    end
    handles.data2.sylEdges = round(handles.data2.sylEdges*handles.data.globalStretch);
end

%Create a linear path between the two spectrograms
ts1 = [1;handles.data1.sylEdges;size(handles.data1.un_meanSpec,2)];
ts2 = [1;handles.data2.sylEdges;size(handles.data2.un_meanSpec,2)];
handles.data.d2dAudioPath = [ts1,ts2];

%Calculate truncation indices for matched comparison
[handles.data1.trncInd,handles.data2.trncInd] = truncSeries(handles.data.d2dAudioPath([1,handles.data.targetInt(1)+1,end],:));

%Calculate harmonic power ratio timeseries
[handles.data1.harmPowRatio,handles.data2.harmPowRatio] = harmonicPower(handles);

%Rewarp day 2 data to match temporal structure
handles.data2.re_meanSpec = alignSeriesSTW(handles.data2.un_meanSpec,handles.data.d2dAudioPath);
if isfield(handles.data2,'neuro')
    handles.data2.re_meanNeuro = alignSeriesSTW(handles.data2.un_meanNeuro,handles.data.d2dAudioPath);
end

%Calculate correlation between traces
a = mat2gray(-10*log10(sum(handles.data1.un_meanSpec(:,handles.data1.trncInd),1)));
b = mat2gray(-10*log10(sum(handles.data2.un_meanSpec(:,handles.data2.trncInd),1)));
c = mat2gray(-10*log10(sum(handles.data2.re_meanSpec(:,handles.data1.trncInd),1)));
handles.data.un_Acorr = winPearson(a,b,str2double(get(handles.edit_corrWinA,'String')));
handles.data.re_Acorr = winPearson(a,c,str2double(get(handles.edit_corrWinA,'String')));
if isfield(handles.data2,'neuro')
    handles.data.un_Ncorr = winPearson(handles.data1.un_meanNeuro(handles.data1.trncInd),handles.data2.un_meanNeuro(handles.data2.trncInd),str2double(get(handles.edit_corrWinN,'String')));
    handles.data.re_Ncorr = winPearson(handles.data1.un_meanNeuro(handles.data1.trncInd),handles.data2.re_meanNeuro(handles.data1.trncInd),str2double(get(handles.edit_corrWinN,'String')));
end

%If present, align unwarped neural traces to each other and calculate the warping path
if isfield(handles.data2,'neuro')
    %Get the scaling values and scale the neural timeseries
    handles.data.maxmax = max(max(handles.data1.un_meanNeuro./mean(handles.data1.un_meanNeuro)),max(handles.data2.un_meanNeuro./mean(handles.data2.un_meanNeuro)));
    handles.data.minmin = min(min(handles.data1.un_meanNeuro./mean(handles.data1.un_meanNeuro)),min(handles.data2.un_meanNeuro./mean(handles.data2.un_meanNeuro)));
    handles.data1.SCun_meanNeuro = mat2gray(handles.data1.un_meanNeuro./mean(handles.data1.un_meanNeuro),[handles.data.minmin,handles.data.maxmax]);
    handles.data2.SCun_meanNeuro = mat2gray(handles.data2.un_meanNeuro./mean(handles.data2.un_meanNeuro),[handles.data.minmin,handles.data.maxmax]);
    
    %Warp scaled timeseries to each other
    [~, p_neuroTS] = DTWFinal(handles.data1.SCun_meanNeuro,handles.data2.SCun_meanNeuro);
    %[~, p_neuroSpecsT] = DTWWeightedNormBandFastest(SCun_av_neuroTS1,SCun_av_neuroTS2);
    [p_neuroTS handles.data2.reNFULL_meanNeuro] = alignSeriesPath(handles.data2.un_meanNeuro,p_neuroTS);    
    
    %Get linear path between Neural Spectrograms using song syllables as anchors
    minSpace = str2num(get(handles.edit_minSpace,'String')); minGradient = str2num(get(handles.edit_minGradient,'String')); checkDist = str2num(get(handles.edit_checkDist,'String'));
    [~,handles.data1.salPoints] = salientPoints(handles.data1.SCun_meanNeuro,minSpace,minGradient,checkDist);%pick out salient points in power trace
    [handles.data2.salPoints] = p_neuroTS(handles.data1.salPoints); %Find the corresponding location in the day 2 trace
    tbN2 = makeLinPath(handles.data1.salPoints,handles.data2.salPoints,handles.data.d2dAudioPath(2:end-1,1));

    %Assemble the linear path for warping the neural traces
    LLneuroPath = [handles.data.d2dAudioPath(2:end-1,1),tbN2];
    handles.data.d2dNeuroPath = [1,1;LLneuroPath;size(handles.data1.un_meanNeuro,2),size(handles.data2.un_meanNeuro,2)];
    handles.data2.reNLL_meanNeuro = alignSeriesSTW(handles.data2.un_meanNeuro,handles.data.d2dNeuroPath);
end

%Plot all data
handles = plotAlignments(handles);

%Set status message
set(handles.text_message,'String','Alignment complete');

guidata(hObject, handles);

function check_newEdges_Callback(hObject, eventdata, handles)

function edit_edgeThresh_Callback(hObject, eventdata, handles)

function edit_edgeThresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_targetStart_Callback(hObject, eventdata, handles)

function edit_targetStart_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_targetEnd_Callback(hObject, eventdata, handles)

function edit_targetEnd_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Audio Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_harmPow_Callback(hObject, eventdata, handles)

function check_harmDiff_Callback(hObject, eventdata, handles)

function check_pitch_Callback(hObject, eventdata, handles)

function edit_pitchCenter_Callback(hObject, eventdata, handles)

function edit_pitchCenter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_numHarm_Callback(hObject, eventdata, handles)

function edit_numHarm_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function check_matchGlobal_Callback(hObject, eventdata, handles)

function edit_smWindow_Callback(hObject, eventdata, handles)

function edit_smWindow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_exportData_Callback(hObject, eventdata, handles)
%This will just export the data associated with the Ali Paper.

%Parse out the necessary values
vernier.targetStart = str2double(get(handles.edit_targetStart,'String'))+1;
vernier.targetEnd = str2double(get(handles.edit_targetEnd,'String'))+1;
vernier.audioPath = handles.data.d2dAudioPath;
vernier.neuroPath = handles.data.d2dNeuroPath;
vernier.unwarpCorr = handles.data.un_Ncorr;
vernier.warpCorr = handles.data.re_Ncorr;
vernier.timeVect = handles.data1.timeVect(handles.data1.trncInd)-handles.data1.sylEdges(handles.data.targetInt(1));

%Request location to save data
[fname,pathLoc] = uiputfile('C:\Users\Tim\Desktop\Ali Datasets Again\*.mat');

%...and save to file
save([pathLoc fname],'vernier')


function push_exportImg_Callback(hObject, eventdata, handles)
%Calculate the shift indices required to align plots tot he selected intervals
targetTimes = [handles.data1.sylEdges(handles.data.targetInt(1)),handles.data2.sylEdges(handles.data.targetInt(1))];

%Set plotting limits
negLim = -1*max(targetTimes);
posLim = max(handles.data.d2dAudioPath(end,1)-targetTimes(1),handles.data.d2dAudioPath(end,2)-targetTimes(2));

%Create Figure
figure

%Plot aligned Spectrograms
subplot(8,1,1)
cla; hold on
imagesc(-targetTimes(1),-0.5,-handles.data1.un_meanSpec); axis xy
barlims = ones(length(handles.data1.sylEdges),1)*[-0.5,90.5];
line([handles.data1.sylEdges-targetTimes(1),handles.data1.sylEdges-targetTimes(1)]',barlims','Color','k','LineStyle','--','LineWidth',1)
axis tight
xlim([negLim,posLim])
set(gca,'XTick',[],'YTick',[],'TickDir','out','Box','off');
title('Baseline')

subplot(8,1,2)
cla; hold on
imagesc(-targetTimes(2),-0.5,-handles.data2.un_meanSpec); axis xy
line([handles.data2.sylEdges-targetTimes(2),handles.data2.sylEdges-targetTimes(2)]',barlims','Color','k','LineStyle','--','LineWidth',1)
axis tight
xlim([negLim,posLim])
xticks = unique([fliplr(-1*(0:100:-negLim)),0:100:posLim]);
set(gca,'XTick',[],'YTick',[],'TickDir','out','Box','off');
title('Post-CAF')

%Plot aligned audio power envelops
subplot(8,1,3)
cla; hold on
plot(handles.data1.timeVect-targetTimes(1),mat2gray(-10*log10(sum(handles.data1.un_meanSpec,1))))
plot(handles.data2.timeVect-targetTimes(2),mat2gray(-10*log10(sum(handles.data2.un_meanSpec,1))),'r')
plot(handles.data1.timeVect-targetTimes(1),mat2gray(-10*log10(sum(handles.data2.re_meanSpec,1))),'g')
line([0,0],[0,1],'Color',[0 0 0],'LineStyle','--','LineWidth',1.5)
axis tight;
xlim([negLim,posLim])
set(gca,'XTick',[],'YTick',[0.5,1],'TickDir','out','Box','off')
title('Normalized Mean Audio Power')

subplot(8,1,4)
cla; hold on
plot(handles.data1.timeVect(handles.data1.trncInd)-targetTimes(1),handles.data.un_Acorr,'k')
plot(handles.data1.timeVect(handles.data1.trncInd)-targetTimes(1),handles.data.re_Acorr,'g')
line([0,0],[-1,1],'Color',[0 0 0],'LineStyle','--','LineWidth',1.5)
axis tight;
xlim([negLim,posLim])
ylim([-1,1.1])
set(gca,'XTick',[],'YTick',[-1,0,1],'TickDir','out','Box','off')
hold off
title('Running Correlation in Audio Power, 50ms Window')

%Plot aligned pitch/harmonic power ratios
subplot(8,1,5)
cla; hold on
if get(handles.check_harmPow,'Value')
    plot(handles.data1.timeVect-targetTimes(1),handles.data1.harmPowRatio)
    plot(handles.data2.timeVect-targetTimes(2),handles.data2.harmPowRatio,'r')
    line([0,0],[0.8,1.2],'Color',[0 0 0],'LineStyle','--','LineWidth',1.5)
    xlim([negLim,posLim])
    ylim([.8,1.2])
    set(gca,'XTick',[],'YTick',[0.85,1,1.15],'TickDir','out','Box','off')
    title(['Harmonic Power Ratio (@ ' num2str(get(handles.edit_pitchCenter,'String')) 'Hz)'])
elseif get(handles.check_harmDiff,'Value')
    diffplot = (handles.data1.harmPowRatio(handles.data1.trncInd)-handles.data2.harmPowRatio(handles.data2.trncInd)).^2;
    area(handles.data1.timeVect(handles.data1.trncInd)-targetTimes(1),diffplot);
    axis tight
    line([0,0],ylim,'Color',[0 0 0],'LineStyle','--','LineWidth',1.5)
    xlim([negLim,posLim])
    set(gca,'XTick',[],'YTick',ylim,'TickDir','out','Box','off')
    title(['Square Difference in Harmonic Power Ratio (@ ' num2str(get(handles.edit_pitchCenter,'String')) 'Hz)'])
elseif get(handles.check_pitch,'Value')
    
end

%Plot aligned neural power envelops
subplot(8,1,6)
cla;
if isfield(handles.data2,'neuro')
    hold on
    maxmax = max(max(handles.data1.un_meanNeuro./mean(handles.data1.un_meanNeuro)),max(handles.data2.un_meanNeuro./mean(handles.data2.un_meanNeuro)));
    minmin = min(min(handles.data1.un_meanNeuro./mean(handles.data1.un_meanNeuro)),min(handles.data2.un_meanNeuro./mean(handles.data2.un_meanNeuro)));
    %plot(tShift1(Indx1),mat2gray(un_av_neuroTS1(Indx1)./mean(un_av_neuroTS1(Indx1)),[minmin,maxmax]))
    plot(handles.data1.timeVect-targetTimes(1),mat2gray(handles.data1.un_meanNeuro./mean(handles.data1.un_meanNeuro),[minmin,maxmax]))
    plot(handles.data2.timeVect-targetTimes(2),mat2gray(handles.data2.un_meanNeuro./mean(handles.data2.un_meanNeuro),[minmin,maxmax]),'r')
    plot(handles.data1.timeVect-targetTimes(1),mat2gray(handles.data2.re_meanNeuro./mean(handles.data2.re_meanNeuro),[minmin,maxmax]),'g')
    line([0,0],[0,1],'Color',[0 0 0],'LineStyle','--','LineWidth',1.5)
    title('Aligned Neural Power Trace')
end
axis tight;
xlim([negLim,posLim])
set(gca,'XTick',[],'YTick',[0.5,1],'TickDir','out','Box','off')
hold off

%Plot neural power correlations
subplot(8,1,7)
cla;
if isfield(handles.data2,'neuro')
    hold on
    plot(handles.data1.timeVect(handles.data1.trncInd)-targetTimes(1),handles.data.un_Ncorr,'k')
    plot(handles.data1.timeVect(handles.data1.trncInd)-targetTimes(1),handles.data.re_Ncorr,'g')
    line([0,0],[-1,1],'Color',[0 0 0],'LineStyle','--','LineWidth',1.5)
end
axis tight;
xlim([negLim,posLim])
ylim([-1,1])
set(gca,'XTick',xticks,'YTick',[-1,0,1],'TickDir','out','Box','off')
hold off
title('Running Correlation in Neural Power, 50ms Window')

%Plot neural power correlations
subplot(8,1,8)
cla; hold on
if isfield(handles.data2,'neuro')
    line([handles.data1.sylEdges-targetTimes(1),handles.data1.sylEdges-targetTimes(1)]',((barlims./3)+0)','Color','b','LineStyle','--','LineWidth',1.5)
    line([handles.data2.sylEdges-targetTimes(2),handles.data2.sylEdges-targetTimes(2)]',((barlims./3)+(barlims(1,2)./3))','Color','r','LineStyle','--','LineWidth',1.5)
    line([handles.data.d2dNeuroPath(2:end-1,2)-targetTimes(1),handles.data.d2dNeuroPath(2:end-1,2)-targetTimes(1)]',((barlims./3)+(2*barlims(1,2)./3))','Color','g','LineStyle','--','LineWidth',1.5)
end
axis tight;
xlim([negLim,posLim])
set(gca,'XTick',[],'YTick',[0.5,1],'TickDir','out','Box','off')
title('Interval Neural Activity Timeseries Alignments')
hold off

function edit_corrWinN_Callback(hObject, eventdata, handles)

function edit_corrWinN_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_corrWinA_Callback(hObject, eventdata, handles)

function edit_corrWinA_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_minSpace_Callback(hObject, eventdata, handles)

function edit_minSpace_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_minGradient_Callback(hObject, eventdata, handles)

function edit_minGradient_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_checkDist_Callback(hObject, eventdata, handles)

function edit_checkDist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
