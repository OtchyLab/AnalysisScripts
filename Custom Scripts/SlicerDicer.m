function varargout = SlicerDicer(varargin)
% SLICERDICER M-file for SlicerDicer.fig
%      SLICERDICER, by itself, creates a new SLICERDICER or raises the existing
%      singleton*.
%
%      H = SLICERDICER returns the handle to a new SLICERDICER or the handle to
%      the existing singleton*.
%
%      SLICERDICER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SLICERDICER.M with the given input arguments.
%
%      SLICERDICER('Property','Value',...) creates a new SLICERDICER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SlicerDicer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SlicerDicer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SlicerDicer

% Last Modified by GUIDE v2.5 05-Jan-2011 13:27:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SlicerDicer_OpeningFcn, ...
                   'gui_OutputFcn',  @SlicerDicer_OutputFcn, ...
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
end


% --- Executes just before SlicerDicer is made visible.
function SlicerDicer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SlicerDicer (see VARARGIN)


%Load program constants via parameters file (default location)
handles.ParamFileLoc = 'C:\Users\Tim\Desktop\parameters.mat';
%handles.ParamFileLoc = '/Volumes/Data Disk/MATLAB code/Test/parameters.mat';

load(handles.ParamFileLoc);

%Copy all parameter values to the handles structure
handles.param_fields = fieldnames(param);
for i=1:length(handles.param_fields)
    eval(['handles.' char(handles.param_fields(i)) '=param.' char(handles.param_fields(i)) ';']);
end

%Copy parameter values to GUI displays
set(handles.edit_MinSylInterval, 'String', num2str(handles.min_syl_pause));
set(handles.edit_MinSylLength, 'String', num2str(handles.min_syl_length));
set(handles.edit_MaxSylLength, 'String', num2str(handles.max_syl_length));
set(handles.edit_MaxSylBuffer, 'String', num2str(handles.syl_buffer));
set(handles.edit_SylThreshGain, 'String', num2str(handles.detectgain));
set(handles.edit_MinSylperSong, 'String', num2str(handles.minsylpersong));
set(handles.edit_th_winsize, 'String', num2str(handles.th_winsize));
set(handles.edit_th_winstep, 'String', num2str(handles.th_winstep));

handles.filenum = 0;

% Choose default command line output for SlicerDicer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SlicerDicer wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = SlicerDicer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = LoadFile(handles, filenum)
% Function loads new song file (i.e., a WAV) into working memory

%Test the requested file for existence/readability
try
    handles.filename = handles.files(filenum).name;
    [handles.wav]=wavread(handles.files(filenum).name);
catch
    error(['File ' handles.filename ' does not exist or is unreadable...continuing to next file.']);
    %set(handles.edit_FileNum,'String',num2str(handles.filenum));
    %uiwait;
    return;
end

%Successfully loaded file, so update displays
handles.filenum = filenum;
set(handles.text_CurFileName, 'String', handles.filename);
%set(handles.text_TotalFiles, 'String', num2str(length(handles.files)));
set(handles.edit_FileNum,'String',num2str(handles.filenum));
handles.time = linspace(0, (length(handles.wav)-1)/handles.fs, length(handles.wav));
if(handles.endNdx - handles.startNdx > .01)
    width = min(handles.endNdx - handles.startNdx, length(handles.wav));
else
    width = length(handles.wav);
end
handles.startNdx = 1;
handles.endNdx = width;
%handles.fs = fs;

[handles] = getAudioPower(handles);

handles.selectedSyll = -1;
handles.rectHandles = [];
handles.txtHandles = [];
handles.loadsuccess = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = LoadBDT(handles, filenum)
% Function loads new song file (i.e., BDT) and annotation into working memory

%Test the requested file for existence/readability
try
    handles.filename = char(handles.files(filenum));
    %[handles.wav]=wavread(handles.files(filenum));
    [HWChannels, handles.wav, time, startSamp] = daq_readDatafileBence([handles.sourcefolder,'\',handles.filename]);
    WavMax = max(handles.wav); WavMin = min(handles.wav);
    ScaleFact = min(abs(.98/WavMax), abs(.98/WavMin));
    handles.wav = ScaleFact*handles.wav; %Scale the WAV values to max +/-1.0
    if max(handles.wav)>1 || min(handles.wav)<-1
        WARNING = 'Will Clip'
    end
catch
    error(['File ' handles.filename ' does not exist or is unreadable...continuing to next file.']);
    %set(handles.edit_FileNum,'String',num2str(handles.filenum));
    %uiwait;
    return;
end

%Successfully loaded file, so update displays
handles.filenum = filenum;
set(handles.text_CurFileName, 'String', handles.filename);
%set(handles.text_TotalFiles, 'String', num2str(length(handles.files)));
set(handles.edit_FileNum,'String',num2str(handles.filenum));
handles.time = linspace(0, (length(handles.wav)-1)/handles.fs, length(handles.wav));
if(handles.endNdx - handles.startNdx > .01)
    width = min(handles.endNdx - handles.startNdx, length(handles.wav));
else
    width = length(handles.wav);
end
handles.startNdx = 1;
handles.endNdx = width;
%handles.fs = fs;

[handles] = getAudioPower(handles);

handles.selectedSyll = -1;
handles.rectHandles = [];
handles.txtHandles = [];
handles.loadsuccess = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function handles = setSelectedSyllable(handles, newSelectedSyll)
%make sure selection is in bounds
%currAnnot = handles.audioAnnotation.get(handles.filename);
if(newSelectedSyll == -1)
    handles.selectedSyll = -1;
    return;
elseif(newSelectedSyll < 1 || newSelectedSyll > length(handles.segments))
    return
end

if(handles.selectedSyll ~= -1 && handles.selectedSyll <= length(handles.segments))
    %unhighlight the old selected syll
    if(handles.segType(handles.selectedSyll) == -1)
        color = 'red';
    else
        color = 'blue';
    end
    if(size(handles.rectHandles,1) >= handles.selectedSyll)
        for(nRect = 1:size(handles.rectHandles,2))
            if(ishandle(handles.rectHandles(handles.selectedSyll,nRect)))
                set(handles.rectHandles(handles.selectedSyll,nRect),'FaceColor',color,'EdgeColor','none');
                if (nRect==2)
                    set(handles.rectHandles(handles.selectedSyll,nRect),'EdgeColor',color);
                end
            end
        end
    end
end

%set the selected syllable to return.
selectedSyll = newSelectedSyll;
handles.selectedSyll = selectedSyll;

%highlight the new selected syllable
if(size(handles.rectHandles,1) >= selectedSyll)
    for(nRect = 1:size(handles.rectHandles,2))
        if(ishandle(handles.rectHandles(selectedSyll,nRect)))
            set(handles.rectHandles(selectedSyll,nRect),'FaceColor','green','EdgeColor','none');
            if (nRect==2)
                set(handles.rectHandles(selectedSyll,nRect),'EdgeColor','green');
            end
        end
    end
end

%recenter if necessary
% filenum = handles.filenum;
% fs = handles.fs;
% startTime = (handles.startNdx-1)/handles.fs;
% endTime = (handles.endNdx-1)/handles.fs;
% syllEndTime = handles.segments(selectedSyll,2);
% syllStartTime = handles.segments(selectedSyll,1);
% if(syllEndTime > endTime || syllStartTime < startTime)
%     width = handles.endNdx - handles.startNdx;
%     mid = (syllStartTime*fs) + 1;
%     handles.startNdx = round(mid - width/2);
%     handles.endNdx = round(mid + width/2);
%     if(handles.startNdx<1)
%         handles.startNdx = 1;
%         handles.endNdx = width + handles.startNdx;
%     end
%     if(handles.endNdx > length(handles.wav))
%         handles.startNdx = length(handles.wav) - width;
%         handles.endNdx = length(handles.wav);
%     end
%     %handles = updateZoomOrPan(handles);
%     handles = drawZoomSyllableRects(handles);
% end
handles = updateZoomOrPan(handles);
handles = drawZoomSyllableRects(handles);
handles = UpdateFeatures(handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = UpdateFeatures(handles)
% Function Executes Feature Extraction from the currently selected syllable
% and displays it in a GUI table

%If a valid Syllable is selected, extract features from that section of wav
if handles.selectedSyll > 0 && handles.selectedSyll <= size(handles.segments,1)
    [handles, FileFeatureSet, FileKeys, FeatureTotals] = Extract33Features(handles, handles.segments(handles.selectedSyll,:));
end

%Format feature data for display on GUI
FileFeatureSet = (round(FileFeatureSet*100))/100; %Truncate values to 2 decimal places
%FileFeatureSet = num2str(FileFeatureSet); %Convert to string
FeatureTotals = (round(FeatureTotals*100))/100; %Truncate values to 2 decimal places
%FeatureTotals = num2str(FeatureTotals); %Convert to string

%Copy data to GUI structures
set(handles.text_SylDurT, 'String', num2str(FileFeatureSet(1,1))); %Duration

set(handles.text_FundFreq1, 'String', num2str(FileFeatureSet(1,2))); %Fundamental Freq
set(handles.text_FundFreq2, 'String', num2str(FileFeatureSet(1,3))); %
set(handles.text_FundFreq3, 'String', num2str(FileFeatureSet(1,4))); %
set(handles.text_FundFreq4, 'String', num2str(FileFeatureSet(1,5))); %
set(handles.text_FundFreqT, 'String', num2str(FeatureTotals(1,1))); %

set(handles.text_FreqSlope1, 'String', num2str(FileFeatureSet(1,7))); %Frequency Slope
set(handles.text_FreqSlope2, 'String', num2str(FileFeatureSet(1,8))); %
set(handles.text_FreqSlope3, 'String', num2str(FileFeatureSet(1,9))); %
set(handles.text_FreqSlope4, 'String', num2str(FileFeatureSet(1,10))); %
set(handles.text_FreqSlopeT, 'String', num2str(FeatureTotals(1,2))); %

set(handles.text_FreqMod1, 'String', num2str(FileFeatureSet(1,30))); %Frequency Mod
set(handles.text_FreqMod2, 'String', num2str(FileFeatureSet(1,31))); %
set(handles.text_FreqMod3, 'String', num2str(FileFeatureSet(1,32))); %
set(handles.text_FreqMod4, 'String', num2str(FileFeatureSet(1,33))); %
set(handles.text_FreqModT, 'String', num2str(FeatureTotals(1,3))); %

set(handles.text_HPAmpT, 'String', num2str(FileFeatureSet(1,6))); %Time to 1/2 Peak

set(handles.text_AmpSlope1, 'String', num2str(FileFeatureSet(1,11))); %Amplitude Slope
set(handles.text_AmpSlope2, 'String', num2str(FileFeatureSet(1,12))); %
set(handles.text_AmpSlope3, 'String', num2str(FileFeatureSet(1,13))); %
set(handles.text_AmpSlopeT, 'String', num2str(FeatureTotals(1,5))); %

set(handles.text_AmpMod1, 'String', num2str(FileFeatureSet(1,26))); %Amplitude Mod
set(handles.text_AmpMod2, 'String', num2str(FileFeatureSet(1,27))); %
set(handles.text_AmpMod3, 'String', num2str(FileFeatureSet(1,28))); %
set(handles.text_AmpMod4, 'String', num2str(FileFeatureSet(1,29))); %
set(handles.text_AmpModT, 'String', num2str(FeatureTotals(1,6))); %

set(handles.text_SEntropy1, 'String', num2str(FileFeatureSet(1,14))); %Spectral Ent
set(handles.text_SEntropy2, 'String', num2str(FileFeatureSet(1,15))); %
set(handles.text_SEntropy3, 'String', num2str(FileFeatureSet(1,16))); %
set(handles.text_SEntropy4, 'String', num2str(FileFeatureSet(1,17))); %
set(handles.text_SEntropyT, 'String', num2str(FeatureTotals(1,7))); %

set(handles.text_TEntropy1, 'String', num2str(FileFeatureSet(1,18))); %Temporal Ent
set(handles.text_TEntropy2, 'String', num2str(FileFeatureSet(1,19))); %
set(handles.text_TEntropy3, 'String', num2str(FileFeatureSet(1,20))); %
set(handles.text_TEntropy4, 'String', num2str(FileFeatureSet(1,21))); %
set(handles.text_TEntropyT, 'String', num2str(FeatureTotals(1,8))); %

set(handles.text_STEntropy1, 'String', num2str(FileFeatureSet(1,22))); %Spectral Temp Ent
set(handles.text_STEntropy2, 'String', num2str(FileFeatureSet(1,23))); %
set(handles.text_STEntropy3, 'String', num2str(FileFeatureSet(1,24))); %
set(handles.text_STEntropy4, 'String', num2str(FileFeatureSet(1,25))); %
set(handles.text_STEntropyT, 'String', num2str(FeatureTotals(1,9))); %

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = UpdateAxes(handles)
% Function updates GUI axes display (does not execute during batch
% processing)

%Plot waveform on the axesNavBar
axes(handles.axesNavBar);
plotTimeSeriesQuick(handles.time, handles.wav);
set(handles.axesNavBar, 'YLim', [-.25, 25]);

%Update zoom/pan position for a;; axes if necessary
%handles = updateZoomOrPan(handles);

%Check to see if file should be segmented on opening; if so, do it
if get(handles.check_SegOnOpen,'Value')
    handles.thresh_done = 0;
    %This is the decision point for using WavChop or using the annotation
    %file for segmentation
    if ~get(handles.check_UseAnnotation,'Value')
        handles = ChopWave(handles);
    else
        %Copy segmentation from annotation file to structure and scale
        handles.segments = handles.fs*[handles.annotation{handles.filenum}.segFileStartTimes', handles.annotation{handles.filenum}.segFileEndTimes'];
        %Create structure for tracking syllable types for the current file
        if ~isempty(handles.segments)
            handles.segType = handles.annotation{handles.filenum}.segType;
        end
    end
    
    %Update zoom/pan position for a;; axes if necessary
    handles = updateZoomOrPan(handles);
   
    %thresholdVector = handles.low_amplitude_th*ones(length(handles.startNdx:handles.endNdx),1);
    
     %axes(handles.axesPower)
     %hold on
     %plotTimeSeriesQuick(handles.startNdx:handles.endNdx, thresholdVector);
     %hold off
    
%     if ~isempty(handles.segments)
%         %Draw syllable boundaries on the display bars
%         handles.rectHandles = [];
%         handles.txtHandles = [];
%         handles = drawNavBarSyllableRects(handles);
%         handles = drawZoomSyllableRects(handles);
%    end
else
    %Update zoom/pan position for a;; axes if necessary
    handles = updateZoomOrPan(handles);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = updateZoomOrPan(handles)
%Function updates the Spectrogram and Power axes

startTime = (handles.startNdx-1)/handles.fs;
endTime = (handles.endNdx-1)/handles.fs;
startNdx = handles.startNdx;
endNdx = handles.endNdx;

%fix navigator rect
axes(handles.axesNavBar);

% if(ishandle(handles.navRect))
%     delete(handles.navRect);
% end

x = [startTime, startTime, endTime, endTime, startTime];
y = ylim;
y = [y(1), y(2), y(2), y(1), y(1)];
handles.navRect = line(x,y);
set(handles.navRect,'Color', 'red');
set(handles.navRect,'LineWidth', 2);
set(handles.navRect,'HitTest','off');

%Replot Sprectrogram and Power
axes(handles.axesSpectrogram);
displaySpecgramQuick(handles.wav(startNdx:endNdx), handles.fs, [0,15000],[],startTime);
xlabel('');
ylabel('');
set(handles.axesSpectrogram,'XTick',[]);

axes(handles.axesPower);
plotTimeSeriesQuick(handles.time(startNdx:endNdx), handles.power(startNdx:endNdx));
set(handles.axesPower,'ButtonDownFcn','');
set(handles.axesPower,'XTick',[]);
%if handles.thresh_done == 1; %isfield(handles, 'thresSyll') && isfield(handles, 'thresTrig') && 
if isfield(handles, 'thresh_done') && handles.thresh_done == 1
    E_thresholdVector = handles.thresSyll(startNdx:endNdx);
    S_thresholdVector = handles.thresTrig(startNdx:endNdx);
    hold on
    Eline = plotTimeSeriesQuick(handles.time(startNdx:endNdx), E_thresholdVector);
    set(Eline, 'Color', 'r');
    hold on
    Sline = plotTimeSeriesQuick(handles.time(startNdx:endNdx), S_thresholdVector);
    set(Sline, 'Color', 'g');
end
hold off

% handles.thresSyll = [];
% handles.thresTrig = [];

linkaxes([handles.axesPower, handles.axesSpectrogram],'x');
set(get(handles.axesNavBar,'Parent'), 'KeyPressFcn', @cb_keypress);
set(handles.axesNavBar, 'ButtonDownFcn', @cb_navbar_click);
set(handles.axesSpectrogram, 'ButtonDownFcn', @cb_spectrogram_click);

handles.rectHandles(:,2:4) = -1;
handles.txtHandles(:,1:3) = -1;
%handles = drawZoomSyllableRects(handles);   %Needs to be re-enabled
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = drawNavBarSyllableRects(handles)

startTime = (handles.startNdx-1)/handles.fs;
endTime = (handles.endNdx-1)/handles.fs;
startNdx = handles.startNdx;
endNdx = handles.endNdx;

    syllStartTimes = (handles.segments(:,1)/handles.fs)';
    syllEndTimes = (handles.segments(:,2)/handles.fs)';
    if isempty(handles.segType)
        syllType = -1*ones(size(handles.segments,1),1);
    else
        syllType = handles.segType;
    end
    x = [syllStartTimes; syllStartTimes; syllEndTimes; syllEndTimes; syllStartTimes];   
    color = zeros(3,length(syllStartTimes));
    color(3,syllType~=-1) = 1;
    color(1,syllType==-1) = 1;
     if(handles.selectedSyll~=-1)
         color(:,handles.selectedSyll) = [0,1,0];
     end
    if ~isempty(syllStartTimes) && ~isempty(syllEndTimes)
    %draw nav bar rectangles.
    axes(handles.axesNavBar)
    lims = ylim;
    y = [lims(1); lims(2); lims(2); lims(1); lims(1)];
        for(nSyll = 1:length(syllStartTimes))
            handles.rectHandles(nSyll,1) = patch(x(:,nSyll),y, color(:,nSyll)', 'EdgeColor', 'none', 'FaceAlpha', .2);
        end
    set(handles.rectHandles(:,1), 'HitTest', 'off');
%    else
%        error('No Syllables Detected in Song File')
%        return
    end
%end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = drawZoomSyllableRects(handles)
%Rescales Spectrogram and Power Axes to fit the selected portions
startTime = (handles.startNdx-1)/handles.fs;
endTime = (handles.endNdx-1)/handles.fs;
startNdx = handles.startNdx;
endNdx = handles.endNdx;

%if(handles.audioAnnotation.containsKey(handles.filename))
%    currAnnot = handles.audioAnnotation.get(handles.filename);
    syllStartTimes = (handles.segments(:,1)/handles.fs)';
    syllEndTimes = (handles.segments(:,2)/handles.fs)';
    midSyll = (syllStartTimes + syllEndTimes) / 2;
    if isempty(handles.segType)
        syllType = -1*ones(size(handles.segments,1),1);
    else
        syllType = handles.segType;
    end
    x = [syllStartTimes; syllStartTimes; syllEndTimes; syllEndTimes; syllStartTimes];   
    color = zeros(3,length(syllStartTimes));
    color(3,syllType~=-1) = 1;
    color(1,syllType==-1) = 1;
    if(handles.selectedSyll~=-1)
        color(:,handles.selectedSyll) = [0,1,0];
    end
 
    %Only the visible rectangles need be drawn in the other axes.
    ndx = find(syllStartTimes < endTime & syllEndTimes> startTime);

    axes(handles.axesSpectrogram)
    lims = ylim;
    y = [lims(1); lims(2); lims(2); lims(1); lims(1)];
    for(nSyll = ndx)
        handles.rectHandles(nSyll, 2) = patch(x(:,nSyll),y, color(:,nSyll)', 'EdgeColor', color(:,nSyll)', 'FaceAlpha', 0);
        if(syllType(nSyll) ~= -1)
            handles.txtHandles(nSyll, 1) = text(midSyll(nSyll),mean(ylim),num2str(syllType(nSyll)));
            set(handles.txtHandles(nSyll,1), 'Color', 'black');
            set(handles.txtHandles(nSyll,1), 'FontSize', 14);
            set(handles.txtHandles(nSyll,1), 'FontWeight', 'bold');
            set(handles.txtHandles(nSyll,1), 'HitTest', 'off'); 
            set(handles.txtHandles(nSyll,1), 'HorizontalAlignment', 'center'); 
        else
            handles.txtHandles(nSyll,1) = -1;
        end
    end

    axes(handles.axesPower)
    lims = ylim;
    y = [lims(1); lims(2); lims(2); lims(1); lims(1)];
    for(nSyll = ndx)
        handles.rectHandles(nSyll, 3) = patch(x(:,nSyll),y, color(:,nSyll)', 'EdgeColor', 'none', 'FaceAlpha', .2);
        if(syllType(nSyll) ~= -1)
            handles.txtHandles(nSyll, 2) = text(midSyll(nSyll),mean(ylim),num2str(syllType(nSyll)));
            set(handles.txtHandles(nSyll,2), 'Color', 'black');
            set(handles.txtHandles(nSyll,2), 'FontSize', 14);
            set(handles.txtHandles(nSyll,2), 'FontWeight', 'bold');
            set(handles.txtHandles(nSyll,2), 'HitTest', 'off'); 
            set(handles.txtHandles(nSyll,2), 'HorizontalAlignment', 'center'); 
        else
            handles.txtHandles(nSyll,2) = -1;
        end
    end  
    
    set(handles.rectHandles(ndx,2:3), 'HitTest', 'off'); 
    handles.rectHandles(handles.rectHandles==0) = -1;
    handles.rectHandles(handles.txtHandles==0) = -1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function cb_navbar_click(hObject, evnt)
handles = guidata(hObject);
fs=handles.fs;
mouseMode = get(get(hObject,'Parent'), 'SelectionType');
clickLocation = get(handles.axesNavBar, 'CurrentPoint');
axes(handles.axesNavBar);
if(strcmp(mouseMode, 'open'))
    %if double click then zoom out
    %disp('open');
    %handles.startNdx = 1;
    %handles.endNdx = length(handles.audio);
    %handles = updateZoomOrPan(handles);
elseif(strcmp(mouseMode, 'alt'))
    %pan length of drag.
elseif(strcmp(mouseMode, 'normal'))
    %zoom 
    rect = rbbox;
    endPoint = get(gca,'CurrentPoint'); 
    point1 = clickLocation(1,1:2);              % extract x and y
    point2 = endPoint(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    offset = abs(point1-point2);         % and dimensions
    l = xlim;
    if((offset(1) / l(2))< .005)
        %if didn't drag a rectangle...
        %recenter on the click location.
        width = handles.endNdx - handles.startNdx;
        mid = (p1(1)*fs) + 1;
        handles.startNdx = round(mid - width/2);
        handles.endNdx = round(mid + width/2);
        if(handles.startNdx<1)
            handles.startNdx = 1;
            handles.endNdx = width + handles.startNdx;
        end
        if(handles.endNdx > length(handles.wav))
            handles.startNdx = length(handles.wav) - width;
            handles.endNdx = length(handles.wav);
        end
        plotTimeSeriesQuick(handles.time, handles.wav);
        %handles = drawNavBarSyllableRects(handles);
        handles = updateZoomOrPan(handles);
        handles = drawZoomSyllableRects(handles);%%%%%%%%
    else
        handles.startNdx = max(1,round((p1(1)*fs) + 1));
        handles.endNdx = min(length(handles.wav), round(((p1(1) + offset(1))*fs) + 1));
        plotTimeSeriesQuick(handles.time, handles.wav);
        %handles = drawNavBarSyllableRects(handles);
        handles = updateZoomOrPan(handles);
        handles = drawZoomSyllableRects(handles);
    end
end
guidata(hObject, handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_spectrogram_click(hObject, evnt)
handles = guidata(hObject);
fs = handles.fs;
filenum = handles.filenum;
mouseMode = get(get(hObject,'Parent'), 'SelectionType');
clickLocation = get(handles.axesSpectrogram, 'CurrentPoint');
axes(handles.axesSpectrogram);

if ~isempty(handles.segments)
    nSelect = find(clickLocation(1,1) >= handles.segments(:,1)./handles.fs & ...
                   clickLocation(1,1) <= handles.segments(:,2)/handles.fs);
    if(length(nSelect) == 1)
        if(nSelect ~= handles.selectedSyll)
            handles = setSelectedSyllable(handles, nSelect);
            
             axes(handles.axesNavBar)
             plotTimeSeriesQuick(handles.time, handles.wav);
             %handles = drawNavBarSyllableRects(handles);
        end

    end
end

guidata(hObject, handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles] = getAudioPower(handles)
%Function reads in WAV file and returns audio power function

%Copy passed in value to local variable
audio = handles.wav;

%De-mean the audio
audio = audio - mean(audio);

%Include only very relevant power:

%Low pass below 12000Hz
filt2.order = 50; %sufficient for 44100Hz or lower
filt2.win = hann(filt2.order+1);
filt2.cutoff = 12000; %Hz
filt2.fs = handles.fs;
filt2.lpf = fir1(filt2.order, filt2.cutoff/(filt2.fs/2), 'low', filt2.win);
audio = filtfilt(filt2.lpf, 1, audio);

%High pass above 300Hz
filt3.order = 50; %sufficient for 44100Hz or lower
filt3.win = hann(filt3.order+1);
filt3.cutoff = 300; %Hz
filt3.fs = handles.fs;
filt3.hpf = fir1(filt3.order, filt3.cutoff/(filt3.fs/2), 'high', filt3.win);
audio = filtfilt(filt3.hpf, 1, audio);

%compute power
audioPow= audio.^2; 

%smooth the power, lpf:
%  filt.order = 100; 
%  filt.win = hann(filt.order+1);
%  filt.cutoff = 50; %Hz
%  filt.fs = handles.fs;
%  filt.lpf = fir1(filt.order, filt.cutoff/(filt.fs/2), 'low', filt.win);
%  audioPow = filtfilt(filt.lpf, 1, audioPow);

% Construct blurring window.
windowWidth = int16(handles.window);
%halfWidth = 500;
halfWidth = windowWidth / 2;
%gaussFilter = gausswin(1000);
gaussFilter = gausswin(handles.window);
gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.

% Do the filtering (x1)
audioPow = conv(audioPow, gaussFilter);
audioPow = audioPow((str2num(int2str(halfWidth))-1):(end-str2num(int2str(halfWidth))));

% Do the filtering (x2)
audioPow = conv(audioPow, gaussFilter);
audioPow = audioPow((str2num(int2str(halfWidth))-1):(end-str2num(int2str(halfWidth))));

% Do the filtering (x3)
audioPow = conv(audioPow, gaussFilter);
audioPow = audioPow((str2num(int2str(halfWidth))-1):(end-str2num(int2str(halfWidth))));

% Do the filtering (x4)
audioPow = conv(audioPow, gaussFilter);
audioPow = audioPow((str2num(int2str(halfWidth))-1):(end-str2num(int2str(halfWidth))));

%Compute log power
audioLogPow = log(audioPow + eps);

%Copy to output structure
handles.wav = audio;
handles.power = audioLogPow;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = ChopWave_old(handles)
%Function segments current file into segments.

%Determine Dynamic Threshold at 50% above mean silence value
handles.segments = [];
[n(1,:),n(2,:)]=hist(handles.power,100);  %Bin m_amplitude data
Max_bin=max(n(1,1:10));
mean_silence=mean(n(2,n(1,:)==Max_bin));  %Find the mode of silent amplitude
handles.low_amplitude_th=mean_silence-(mean_silence*handles.detectgain/100); %Set threashold this % above silence value

%Threshold to detect syllables
a = [0; handles.power; 0];
handles.segments(:,1) = find(a(1:end-1)<handles.low_amplitude_th & a(2:end)>=handles.low_amplitude_th);
handles.segments(:,2) = find(a(1:end-1)>=handles.low_amplitude_th & a(2:end)<handles.low_amplitude_th);
a = a(2:end-1);

% Eliminate short intervals between 'detected' syllables
j = [find(handles.segments(2:end,1)-handles.segments(1:end-1,2) < handles.min_syl_pause*handles.fs/1000); length(handles.segments)];
%j = [find(handles.segments(:,1)-handles.segments(:,2) < handles.min_syl_pause*handles.fs/1000); length(handles.segments)];
m=0;
while length(j)>1 && m<75
    if j(1)~=1 && j(1)~=length(handles.segments)
        if handles.segments(j(1),1)-handles.segments(j(1)-1,2)<= handles.segments(j(1),2)-handles.segments(j(1)+1,1)
            handles.segments(j(1)-1,2)=handles.segments(j(1),2); %Merge with previous syllable
            handles.segments = [handles.segments(1:j(1)-1,:); handles.segments(j(1)+1:end,:)];
        else
            handles.segments(j(1)+1,1)=handles.segments(j(1),1); %Merge with next syllable
            handles.segments = [handles.segments(1:j(1)-1,:); handles.segments(j(1)+1:end,:)];
        end
    end
j = [find(handles.segments(2:end,1)-handles.segments(1:end-1,2) < handles.min_syl_pause*handles.fs/1000); length(handles.segments)];
m=m+1;
end

i = [find(handles.segments(2:end,1)-handles.segments(1:end-1,2) > handles.min_syl_pause*handles.fs/1000); length(handles.segments)];
%i = [find(handles.segments(:,1)-handles.segments(:,2) > handles.min_syl_pause*handles.fs/1000); length(handles.segments)];

try
    handles.segments = [handles.segments([1; i(1:end-1)+1],1) handles.segments(i,2)];
catch me
    i = 1:(length(handles.segments)-1);
    handles.segments = [handles.segments([1; i(1:end-1)+1],1) handles.segments(i,2)];
end

% Eliminate short syllables
%i = find(handles.segments(:,2)-handles.segments(:,1) > handles.min_syl_length*handles.fs/1000);
i = find(handles.segments(2:end,2)-handles.segments(1:end-1,1) > handles.min_syl_length*handles.fs/1000);
%handles.segments = handles.segments(i,:);
handles.segments = [handles.segments(i,1) handles.segments(i+1,2)];

% Eliminate any detected syllable that coincides with the end of the wav file
if handles.segments(end,2) >= (length(handles.power)-1000)
    handles.segments = handles.segments(1:end-1,:);
end

% Attach dynamic buffer
for i=1:size(handles.segments,1)
    %Attach buffer to start position
    if handles.segments(i,1)>round(handles.syl_buffer*handles.fs/1000)
        search = (handles.segments(i,1)-round(handles.syl_buffer*handles.fs/1000)):(handles.segments(i,1));
        [min_val, ind] = min(a(search));
        if i~=1
            handles.segments(i,1) = max(ind+(handles.segments(i,1)-round(handles.syl_buffer*handles.fs/1000))-1, handles.segments(i-1,2));
        else
            handles.segments(i,1) = ind+(handles.segments(i,1)-round(handles.syl_buffer*handles.fs/1000))-1;
        end
    end
    %Attach buffer to end position
    if handles.segments(i,2)<(length(a)-round(handles.syl_buffer*handles.fs/1000))
        search = (handles.segments(i,2):(handles.segments(i,2)+round(handles.syl_buffer*handles.fs/1000)));
        [min_val, ind] = min(a(search));
        if i~=size(handles.segments,1)
            handles.segments(i,2) = min(ind+(handles.segments(i,2))+1, handles.segments(i+1,1));
        else
            handles.segments(i,2) = ind+(handles.segments(i,2)+1);
        end
    end
end

%Create structure for tracking syllable types for the current file
%Initialize with -1 for "unsorted"
handles.segType = -1*ones(size(handles.segments,1),1);

end %function end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = ChopWave(handles)
%This is a second method for cutting a complete songfile in to individual
%syllables.  The function utilizes several thresholds and criteria:

%1. syllable trig threshold power - the syllable has to reach this
%   threshold at least once.
%2. syllable cont threshold power - all time between this power is
%   considered the syllable.
%3. min interval duration.
%4. min syllable duration.
%5. max syllable duration

%audio is the rawaudio file
%fs is the sampling rate.
%syllStartTimes: the start time of each syllable.  This is not the onset of
%sound, but the time at which the audio level crosses above threshold.
%syllEndTimes: the end time of each syllable.  This is not the offset of
%sound, but the time at which the audio level crosses below threshold.
%noiseEst, noiseStd: the estimated noise level and variance.
%soundEst, soundStd: the estimated signal level and variance.

%Set Chop Flag
handles.thresh_done = 0;

syllStartTimes = [];
syllEndTimes = [];
filesForManual = [];
AudioWindows = [];

fMinSyllDuration = handles.min_syl_length*handles.fs/1000;
fMinIntervalDuration = handles.min_syl_pause*handles.fs/1000;
fMaxSyllDuration = handles.max_syl_length*handles.fs/1000;
threshold = handles.detectgain;
trigger = handles.syl_buffer;
win_size = floor((handles.th_winsize/1000)*handles.fs);
win_step = floor((handles.th_winstep/1000)*handles.fs);

%Window the power vector by the input parameters to dynamically calculate
%the threshold
startPos = 1;
SmallestWindow = floor(win_size/2);

while startPos<length(handles.power)
    if (startPos+win_size)>(length(handles.power)-4410)
        AudioWindows = [AudioWindows; startPos length(handles.power)];
        startPos = length(handles.power) +  1;
    else
        AudioWindows = [AudioWindows; startPos startPos+win_size];
        startPos = startPos + win_step;
    end
end
        
numWindows = size(AudioWindows,1);
thresSyllwins = zeros(numWindows,length(handles.power));
SyllKey = zeros(numWindows,length(handles.power));
thresTrigwins = zeros(numWindows,length(handles.power));
TrigKey = zeros(numWindows,length(handles.power));

for i = 1:numWindows
    %Estimate threshold to discrimate between sound and noise using a mixture
    %of two gaussians model.
    [noiseEst, soundEst, noiseStd, soundStd] = aSAP_estimateTwoMeans(handles.power(AudioWindows(i,1):AudioWindows(i,2)));    
    if(noiseEst>soundEst)
        return
    end
        %Compute the optimal classifier between the two gaussians...
        p(1) = 1/(2*soundStd^2) - 1/(2*noiseStd^2);
        p(2) = (noiseEst)/(noiseStd^2) - (soundEst)/(soundStd^2);
        p(3) = (soundEst^2)/(2*soundStd^2) - (noiseEst^2)/(2*noiseStd^2) + log(soundStd/noiseStd);
        disc = roots(p);
        disc = disc(disc>noiseEst & disc<soundEst);
       if(isempty(disc))
         return;
       end
        disc = disc(1);

    %Set the thresholds based on these estimates
    thresSyllwins(i,AudioWindows(i,1):AudioWindows(i,2)) = noiseEst + threshold * (disc - noiseEst); %threshold for the edge of a syllable
    Key(i,AudioWindows(i,1):AudioWindows(i,2)) = 1;
    thresTrigwins(i,AudioWindows(i,1):AudioWindows(i,2)) = soundEst - trigger * (soundEst - disc); % it is only a syllable if the power is above this value (get rid of noise)
end



if size(thresSyllwins,1)>1
    handles.thresTrig = smooth(sum(thresSyllwins)./sum(Key),4410); %Smooth with a 100ms window
    handles.thresSyll = smooth(sum(thresTrigwins)./sum(Key),4410);
else
    handles.thresTrig = smooth(thresSyllwins./Key,4410); %Smooth with a 100ms window
    handles.thresSyll = smooth(thresTrigwins./Key,4410);
end

%Find threshold crossings:
% [trigCross, junk] = detectThresholdCrossings(handles.power, handles.thresTrig, true);
% [syllUpCross, syllDownCross] = detectThresholdCrossings(handles.power, handles.thresSyll, true);

%Must modify the below routines to handle the threshold vectors rather than
%threshold scalars
[trigCross, junk] = detectThresholdVectCrossings(handles.power, handles.thresTrig, true);
[syllUpCross, syllDownCross] = detectThresholdVectCrossings(handles.power, handles.thresSyll, true);


if(~isempty(syllUpCross) || ~isempty(syllDownCross))
    
    %Eliminated extraneous end crossing...
    if(syllUpCross(1) == 1)
        syllUpCross = syllUpCross(2:end);
    end
    if(syllDownCross(end) == length(handles.power))
        syllDownCross = syllDownCross(1:end-1);
    end

    %Determine syllables present
    nSyll = 0;
    beginSyll = [];
    endSyll = [];
    for(nTrig = 1:length(trigCross))
        up = find(syllUpCross<trigCross(nTrig));
        down = find(syllDownCross>trigCross(nTrig));
        if((length(up)>0) & (length(down)>0))
            nSyll = nSyll + 1;
            beginSyll(nSyll) = syllUpCross(up(end));
            endSyll(nSyll) = syllDownCross(down(1));
        end
    end

    if(length(beginSyll) > 2)
        %Remove small intervals
        intervals = (beginSyll(2:end) - endSyll(1:end-1));
        realGapNdx = find(intervals > fMinIntervalDuration);
        beginSyll = beginSyll([1,realGapNdx+1]);
        endSyll = endSyll([realGapNdx,length(endSyll)]);
    end
        
    %Remove syllables that are too short or long
    durations = (endSyll - beginSyll);
    realSyll = find((durations > fMinSyllDuration) & (durations < fMaxSyllDuration));
    beginSyll = beginSyll(realSyll);
    endSyll = endSyll(realSyll);

    syllStartTimes = (beginSyll -1);
    syllEndTimes = (endSyll-1);
end

%Copy to main segmentation output structure
handles.segments = [syllStartTimes', syllEndTimes'];
if ~isempty(handles.segments)
    % Eliminate any detected syllable that coincides with the end of the wav file
    if handles.segments(end,2) >= (length(handles.power)-1000)
        handles.segments = handles.segments(1:end-1,:);
    end

    %Create structure for tracking syllable types for the current file
    %Initialize with -1 for "unsorted"
    handles.segType = -1*ones(size(handles.segments,1),1);
end
%Set Chop Flag
handles.thresh_done = 1;

end

function [leadingEdgeNdx, fallingEdgeNdx] = detectThresholdVectCrossings(sig, fThres, bAbove)
%Returns indices of values above or below threshold.  LeadingEdge indices
%are those that first cross the thres.  FallingEdge contains the indices of
%the last value to be above the thres.

leadingEdgeNdx = [];
fallingEdgeNdx = [];
if(bAbove)
    exceedsThres = find(sig>fThres);
else
    exceedsThres = find(sig<fThres);
end
    
if ~isempty(exceedsThres)
    ndx = find(diff(exceedsThres)>1);
    leadingEdgeNdx = exceedsThres(ndx + 1);
    leadingEdgeNdx = [exceedsThres(1); leadingEdgeNdx];
    fallingEdgeNdx = exceedsThres(ndx);  
    fallingEdgeNdx = [fallingEdgeNdx; exceedsThres(end)];
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = ExportFileSyls(handles)
%Function saves the detected syllables in the file to a folder

%Check to see if Export folder is valid
if ~isdir(handles.SylExpFolder)
    error('No valid export location selected')
    return
end

%Check to see if there are any syllables detected
if isempty(handles.segments)
    error('No syllables detected in current file')
    return
else  %...otherwise, export syllables as individual wavs
    for c = 1:size(handles.segments,1)
        fname = regexprep(handles.filename, '.wav', ''); %remove lowercase suffix
        fname = regexprep(fname, '.WAV', ''); %remove uppercase suffix
        fname = regexprep(fname, '.BDT', ''); %remove uppercase suffix
        if get(handles.check_UseAnnotation,'Value')
            fname = [fname '_syll_' num2str(c,'%03.f') '_' num2str(handles.annotation{handles.filenum}.segType(c)) '.wav']; %name new file with syl # attached 
            if ~isdir([handles.SylExpFolder '\' num2str(handles.annotation{handles.filenum}.segType(c))])
                mkdir([handles.SylExpFolder '\' num2str(handles.annotation{handles.filenum}.segType(c))])
            end
                wavwrite(handles.wav(round(handles.segments(c,1):min([handles.segments(c,2) length(handles.wav)]))),handles.fs,16,[handles.SylExpFolder '\' num2str(handles.annotation{handles.filenum}.segType(c)) '\' fname]);
        else
            fname = [fname '_syll_' num2str(c,'%03.f') '.wav']; %name new file with syl # attached
            wavwrite(handles.wav(round(handles.segments(c,1):min([handles.segments(c,2) length(handles.wav)]))),handles.fs,16,[handles.SylExpFolder '\' fname]);
        end
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GUI Object function calls      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_FileNum_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FileNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FileNum as text
%        str2double(get(hObject,'String')) returns contents of edit_FileNum as a double
filenum = get(handles.edit_FileNum, 'String');
filenum = str2num(filenum);
handles.startNdx = -1;
handles.endNdx = -1;

%if handles.dataType == 0
    handles = LoadFile(handles, filenum);
%elseif handles.dataType == 1
%    handles = LoadBDT(handles, filenum);
%end

handles = UpdateAxes(handles);
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function edit_FileNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FileNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% --- Executes on button press in push_PrevFile.
function push_PrevFile_Callback(hObject, eventdata, handles)
% hObject    handle to push_PrevFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filenum = handles.filenum - 1;
handles.startNdx = -1;
handles.endNdx = -1;

if handles.dataType == 0
    handles = LoadFile(handles, filenum);
elseif handles.dataType == 1
    handles = LoadBDT(handles, filenum);
end

handles = UpdateAxes(handles);

guidata(hObject, handles);
end

% --- Executes on button press in push_NextFile.
function push_NextFile_Callback(hObject, eventdata, handles)
% hObject    handle to push_NextFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filenum = handles.filenum + 1;
handles.startNdx = -1;
handles.endNdx = -1;

%if handles.dataType == 0
handles = LoadFile(handles, filenum);
%elseif handles.dataType == 1
%    handles = LoadBDT(handles, filenum);
%end

handles = UpdateAxes(handles);

guidata(hObject, handles);
end

% --- Executes on button press in push_SaveParams.
function push_SaveParams_Callback(hObject, eventdata, handles)
% hObject    handle to push_SaveParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Copy all parameter values to the handles structure
for i=1:length(handles.param_fields)
    eval(['param.' char(handles.param_fields(i)) '=handles.' char(handles.param_fields(i))]);
end

if get(handles.check_NewParamFile, 'Value')==1
    handles.ParamFileSaveLoc = uigetdir(pwd,'Parameter Save Location');
    handles.ParamFileSaveLoc = [handles.ParamFileSaveLoc '/Parameters.mat'];
else
    handles.ParamFileSaveLoc = handles.ParamFileLoc;
end

save(handles.ParamFileSaveLoc, 'param')

guidata(hObject, handles);
end

% --- Executes on button press in check_NewParamFile.
function check_NewParamFile_Callback(hObject, eventdata, handles)
% hObject    handle to check_NewParamFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_NewParamFile
end

% --- Executes on button press in push_LoadParams.
function push_LoadParams_Callback(hObject, eventdata, handles)
% hObject    handle to push_LoadParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get user input for which paramters file to import
%NOTE: at this point, it will overwrite the entire handles structure
handles.ParamFileLoc = uigetfile('*.mat', 'Parameters File');
load(handles.ParamFileLoc)

%Copy all parameter values to the handles structure
handles.param_fields = fieldnames(param);
for i=1:length(handles.param_fields)
    eval(['handles.' char(handles.param_fields(i)) '=param.' char(handles.param_fields(i))]);
end

set(handles.edit_MinSylInterval, 'String', num2str(handles.min_syl_pause));
set(handles.edit_MinSylLength, 'String', num2str(handles.min_syl_length));
set(handles.edit_MaxSylBuffer, 'String', num2str(handles.syl_buffer));
set(handles.edit_SylThreshGain, 'String', num2str(handles.detectgain));
set(handles.edit_MinSylperSong, 'String', num2str(handles.minsylpersong));

guidata(hObject, handles);
end

% --- Executes on button press in push_OpenFolder.
function push_OpenFolder_Callback(hObject, eventdata, handles)
% hObject    handle to push_OpenFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Request user input for which folder to read WAVs from
handles.sourcefolder = uigetdir(pwd,'Source Data Folder');

%songsource = 'C:\Documents and Settings\Bird Brain\Desktop\Test\Wave Chopping Test Folder\Source Data\';
cd(handles.sourcefolder);
handles.files = dir('*.wav');
%handles.files = [handles.files; dir('*.wav')];

if ~isempty(handles.files)
    filenum = 1;
    handles.startNdx = -1;
    handles.endNdx = -1;
    set(handles.text_TotalFiles, 'String', num2str(length(handles.files)));
    handles = LoadFile(handles, filenum);
    handles = UpdateAxes(handles);
else
    error('No WAVs found in the selected folder')
    return
end

%Set flag for loading WAVs (and not BDTs)
handles.dataType = 0;
handles.thresh_done = 0;
guidata(hObject, handles);
end



% --- Executes on button press in push_OpenBDT.
function push_OpenBDT_Callback(hObject, eventdata, handles)
% hObject    handle to push_OpenBDT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Request user input for which folder to read BDTs from
handles.sourcefolder = uigetdir(pwd,'Source Data Folder');

%songsource = 'C:\Documents and Settings\Bird Brain\Desktop\Test\Wave Chopping Test Folder\Source Data\';
cd(handles.sourcefolder);
files = dir('*.BDT');
%handles.files = [handles.files; dir('*.wav')];

if isempty(files)
    error('No BDTs found in the selected folder')
    return
end

%Request user input for annatation file
[File,Path] = uigetfile(pwd,'Annotation File');
% Path = 'C:\Users\Tim\Desktop\Pur238 Sorted Cells\Pur238 063010 107dph\';
% File = 'Pur238_20100630_2765u_annotation.mat';
% handles.sourcefolder = Path;
%Load annotated file list
load([Path, File], 'keys');
handles.files = keys;

%Load annotation elements
load([Path, File], 'elements');
handles.annotation = elements;

if ~isempty(handles.files) || ~isempty(handles.annotation)
    filenum = 1;
    handles.startNdx = -1;
    handles.endNdx = -1;
    set(handles.text_TotalFiles, 'String', num2str(length(handles.annotation)));
    handles = LoadBDT(handles, filenum);
    handles = UpdateAxes(handles);
else
    error('Annotation file failed to load or was not properly structured')
    return
end

%Set flag for loading BDTs (and not WAVs)
handles.dataType = 1;

guidata(hObject, handles);
end

% --- Executes on button press in check_UseAnnotation.
function check_UseAnnotation_Callback(hObject, eventdata, handles)
% hObject    handle to check_UseAnnotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_UseAnnotation
end


% --- Executes on button press in check_SegOnOpen.
function check_SegOnOpen_Callback(hObject, eventdata, handles)
% hObject    handle to check_SegOnOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_SegOnOpen
end

function edit_MinSylInterval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MinSylInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MinSylInterval as text
%        str2double(get(hObject,'String')) returns contents of edit_MinSylInterval as a double

handles.min_syl_pause = get(handles.edit_MinSylInterval, 'String');
handles.min_syl_pause = str2num(handles.min_syl_pause);
startNdx = handles.startNdx;
endNdx = handles.endNdx;
handles = LoadFile(handles, handles.filenum);
handles.startNdx = startNdx;
handles.endNdx= endNdx;
handles = UpdateAxes(handles);

guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function edit_MinSylInterval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MinSylInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit_MinSylLength_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MinSylLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MinSylLength as text
%        str2double(get(hObject,'String')) returns contents of edit_MinSylLength as a double

handles.min_syl_length = get(handles.edit_MinSylLength, 'String');
handles.min_syl_length = str2num(handles.min_syl_length);
startNdx = handles.startNdx;
endNdx = handles.endNdx;
handles = LoadFile(handles, handles.filenum);
handles.startNdx = startNdx;
handles.endNdx= endNdx;
handles = UpdateAxes(handles);

guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function edit_MinSylLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MinSylLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit_MaxSylBuffer_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MaxSylBuffer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MaxSylBuffer as text
%        str2double(get(hObject,'String')) returns contents of edit_MaxSylBuffer as a double

handles.syl_buffer = get(handles.edit_MaxSylBuffer, 'String');
handles.syl_buffer = str2num(handles.syl_buffer);
startNdx = handles.startNdx;
endNdx = handles.endNdx;
handles = LoadFile(handles, handles.filenum);
handles.startNdx = startNdx;
handles.endNdx= endNdx;
handles = UpdateAxes(handles);

guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function edit_MaxSylBuffer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MaxSylBuffer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit_SylThreshGain_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SylThreshGain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SylThreshGain as text
%        str2double(get(hObject,'String')) returns contents of edit_SylThreshGain as a double

handles.detectgain = get(handles.edit_SylThreshGain, 'String');
handles.detectgain = str2num(handles.detectgain);
startNdx = handles.startNdx;
endNdx = handles.endNdx;
handles = LoadFile(handles, handles.filenum);
handles.startNdx = startNdx;
handles.endNdx= endNdx;
handles = UpdateAxes(handles);

guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function edit_SylThreshGain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SylThreshGain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in push_ExportCurSylls.
function push_ExportCurSylls_Callback(hObject, eventdata, handles)
% hObject    handle to push_ExportCurSylls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = ExportFileSyls(handles);
guidata(hObject, handles);

end

% --- Executes on button press in check_ExpAllSyls.
function check_ExpAllSyls_Callback(hObject, eventdata, handles)
% hObject    handle to check_ExpAllSyls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_ExpAllSyls
end

% --- Executes on button press in push_SetExportLoc.
function push_SetExportLoc_Callback(hObject, eventdata, handles)
% hObject    handle to push_SetExportLoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Request user input for where to export wavs to
handles.SylExpFolder = uigetdir(pwd,'Syllable Export Folder');
guidata(hObject, handles);
end

% --- Executes on button press in push_FeatureFileLoc.
function push_FeatureFileLoc_Callback(hObject, eventdata, handles)
% hObject    handle to push_FeatureFileLoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Request user input for where to save Feature File
[handles.FeatureFileName, handles.FeatureFileLoc]=uiputfile('features.mat','Save features as:');

guidata(hObject, handles);

end


% --- Executes on button press in push_ProcessWAVs.
function push_ProcessWAVs_Callback(hObject, eventdata, handles)
% hObject    handle to push_ProcessWAVs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.dataType ~=0
    error('WAV files not loaded.  Are BDTs?')
    return
end

if  get(handles.check_ExpAllSyls, 'Value')==0 && get(handles.check_CreateFeatureFile, 'Value')==0
    error('No processing options chosen')
    return
end

if ~isdir(handles.SylExpFolder) && get(handles.check_ExpAllSyls, 'Value')==1
    error('No valid syllable export location selected')
    return
end

if (~isdir(handles.FeatureFileLoc) || isempty(handles.FeatureFileName)) && get(handles.check_CreateFeatureFile, 'Value')==1
    error('No valid feature file location/name selected')
    return
end

FeatureSet = [];
FeatureKeys = [];
FeatureTotals = [];
keys = [];
m=1;

%Copy all parameter values to the handles structure
for i=1:length(handles.param_fields)
    eval(['settings.' char(handles.param_fields(i)) '=handles.' char(handles.param_fields(i)) ';']);
end

%Process each of the files in the folder in order
for i = 1:length(handles.files)
    if get(handles.check_UseTimeWindow, 'Value')==1
        hour = regexp(handles.files(i).name, '_', 'split'); %Parse file name for hour
        hour=str2double(hour(6));
        if hour>=handles.ProStartTime && hour<=handles.ProEndTime;
            DoIt = 1;
        else
            DoIt = 0;
        end
    else
        DoIt = 1;
    end
    
    if DoIt
        filenum = i;
        handles.startNdx = -1;
        handles.endNdx = -1;
        handles.loadsuccess = 0;
        handles = LoadFile(handles, filenum); %Load file into working memory
        if handles.loadsuccess == 1
            handles = ChopWave(handles); %Segment wave into syllable
            if handles.thresh_done == 1 %The if statement added 5-17-2011 to deal with kick-outs from ChopWave
            if size(handles.segments, 1) >= handles.minsylpersong %See if this file counts as a 'song'
                %Execute if syllables are to be exported
                if get(handles.check_ExpAllSyls, 'Value')==1 
                    handles = ExportFileSyls(handles);
                end
                %Execute if Feature File to be created
                if get(handles.check_CreateFeatureFile, 'Value')==1
                    [handles, FileFeatureSet, FileKeys, FileFeatureTotals] = Extract33Features(handles, handles.segments);
                    FeatureSet = [FeatureSet; FileFeatureSet];
                    FeatureKeys = [FeatureKeys; FileKeys'];
                    FeatureTotals = [FeatureTotals; FileFeatureTotals];
                end
                %Execute if Annotation to be created
                if get(handles.check_CreateAnnotFile, 'Value')==1
                    %Construct keys structure
                    keys(m) = handles.files(i).name;
                    %Construct elements structure
                    elements(m).filenum = m;
                    elements(m).filename = handles.files(i).name;
                    elements(m).segFileStartTimes = handles.segments(:,1);
                    elements(m).segFileEndTimes = handles.segments(:,2);
                    elements(m).segType = handles.segType;
                    elements(m).fs = handles.fs;
                    elements(m).settings = settings;
                    elements(m).sylFeatures = FileFeatureSet;
                    m=m+1;
                end
            end
            end
        end
    end
end

%Save the Syllable features data
if get(handles.check_CreateFeatureFile,'Value')==1
    save([handles.FeatureFileLoc '\' handles.FeatureFileName], 'FeatureSet', 'FeatureKeys', 'FeatureTotals')
end

%Save the Annotation Structures
if get(handles.check_CreateAnnotFile, 'Value')==1
    save([handles.FeatureFileLoc '\' handles.FeatureFileName], 'elements', 'keys')
end

end

% --- Executes on button press in push_ProcessBDTs.
function push_ProcessBDTs_Callback(hObject, eventdata, handles)
% hObject    handle to push_ProcessBDTs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.dataType ~=1
    error('BDT files not loaded.  Are WAVs?')
    return
end

if  get(handles.check_ExpAllSyls, 'Value')==0 && get(handles.check_CreateFeatureFile, 'Value')==0
    error('No processing options chosen')
    return
end

if ~isdir(handles.SylExpFolder) && get(handles.check_ExpAllSyls, 'Value')==1
    error('No valid syllable export location selected')
    return
end

if (~isdir(handles.FeatureFileLoc) || isempty(handles.FeatureFileName)) && get(handles.check_CreateFeatureFile, 'Value')==1
    error('No valid feature file location/name selected')
    return
end

if get(handles.check_UseTimeWindow, 'Value')==1
    error('Cannot (yet) restrict to time window for BDT processing.  Uncheck the box before proceding.')
    return
end


FeatureSet = [];
FeatureKeys = [];
FeatureTotals = [];
SylTypes = [];
keys = [];
m=1;

%Copy all parameter values to the handles structure
for i=1:length(handles.param_fields)
    eval(['settings.' char(handles.param_fields(i)) '=handles.' char(handles.param_fields(i)) ';']);
end

%Process each of the files in the folder in order
for i = 1:length(handles.files)
    if get(handles.check_UseRecWindow, 'Value')
        recnum = regexp(handles.files(i), '_', 'split'); %Parse file name for hour
        recnum=str2double(recnum{1}{2}(2:end));
        if recnum>=handles.ProStartTime && recnum<=handles.ProEndTime;
            DoIt = 1;
        else
            DoIt = 0;
        end
    else
        DoIt = 1;
     end
    
    if DoIt
        filenum = i;
        handles.startNdx = -1;
        handles.endNdx = -1;
        handles.loadsuccess = 0;
        handles = LoadBDT(handles, filenum); %Load file into working memory
        if handles.loadsuccess == 1
            %Copy segmentation from annotation file to structure and scale
            handles.segments = round(handles.fs*[handles.annotation{handles.filenum}.segFileStartTimes', handles.annotation{handles.filenum}.segFileEndTimes']);
            %Create structure for syllable types for the current file
            if ~isempty(handles.segments)
                handles.segType = handles.annotation{handles.filenum}.segType;
            end
            if size(handles.segments, 1) >= handles.minsylpersong %See if this file counts as a 'song'
                
                %Execute if syllables are to be exported
                if get(handles.check_ExpAllSyls, 'Value')==1 
                    handles = ExportFileSyls(handles);
                end
                
                %Execute if Feature File to be created
                if get(handles.check_CreateFeatureFile, 'Value')==1
                    [handles, FileFeatureSet, FileKeys, FileFeatureTotals] = Extract33Features(handles, handles.segments);
                    FeatureSet = [FeatureSet; FileFeatureSet];
                    FeatureKeys = [FeatureKeys; FileKeys'];
                    FeatureTotals = [FeatureTotals; FileFeatureTotals];
                    SylTypes = [SylTypes; handles.segType];
                end
                
                %Execute if Annotation to be created
%                 if get(handles.check_CreateAnnotFile, 'Value')==1
%                     %Construct keys structure
%                     keys(m) = handles.files(i).name;
%                     %Construct elements structure
%                     elements(m).filenum = m;
%                     elements(m).filename = handles.files(i).name;
%                     elements(m).segFileStartTimes = handles.segments(:,1);
%                     elements(m).segFileEndTimes = handles.segments(:,2);
%                     elements(m).segType = handles.segType;
%                     elements(m).fs = handles.fs;
%                     elements(m).settings = settings;
%                     elements(m).sylFeatures = FileFeatureSet;
%                     m=m+1;
%                 end
            end
        end
    end
end

%Save the Syllable features data
if get(handles.check_CreateFeatureFile,'Value')==1
    save([handles.FeatureFileLoc '\' handles.FeatureFileName], 'FeatureSet', 'FeatureKeys', 'FeatureTotals', 'SylTypes')
end

%Save the Annotation Structures
% if get(handles.check_CreateAnnotFile, 'Value')==1
%     save([handles.FeatureFileLoc '\' handles.FeatureFileName], 'elements', 'keys')
% end

end

% --- Executes on button press in check_CreateFeatureFile.
function check_CreateFeatureFile_Callback(hObject, eventdata, handles)
% hObject    handle to check_CreateFeatureFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_CreateFeatureFile
end

function edit_MinSylperSong_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MinSylperSong (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MinSylperSong as text
%        str2double(get(hObject,'String')) returns contents of edit_MinSylperSong as a double

handles.minsylpersong = get(handles.edit_MinSylperSong, 'String');
handles.minsylpersong = str2num(handles.minsylpersong);

guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function edit_MinSylperSong_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MinSylperSong (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in check_CreateAnnotFile.
function check_CreateAnnotFile_Callback(hObject, eventdata, handles)
% hObject    handle to check_CreateAnnotFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_CreateAnnotFile
end

% --- Executes on button press in check_UseTimeWindow.
function check_UseTimeWindow_Callback(hObject, eventdata, handles)
% hObject    handle to check_UseTimeWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_UseTimeWindow
end


function edit_StartTime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_StartTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_StartTime as text
%        str2double(get(hObject,'String')) returns contents of edit_StartTime as a double

handles.ProStartTime = get(handles.edit_StartTime, 'String');
handles.ProStartTime = str2num(handles.ProStartTime);

guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function edit_StartTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_StartTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit_EndTime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_EndTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_EndTime as text
%        str2double(get(hObject,'String')) returns contents of edit_EndTime as a double

handles.ProEndTime = get(handles.edit_EndTime, 'String');
handles.ProEndTime = str2num(handles.ProEndTime);

guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function edit_EndTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_EndTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit_MaxSylLength_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MaxSylLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MaxSylLength as text
%        str2double(get(hObject,'String')) returns contents of edit_MaxSylLength as a double
handles.max_syl_length = get(handles.edit_MaxSylLength, 'String');
handles.max_syl_length = str2num(handles.max_syl_length);
startNdx = handles.startNdx;
endNdx = handles.endNdx;
handles = LoadFile(handles, handles.filenum);
handles.startNdx = startNdx;
handles.endNdx= endNdx;
handles = UpdateAxes(handles);

guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function edit_MaxSylLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MaxSylLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit_th_winsize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_th_winsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_th_winsize as text
%        str2double(get(hObject,'String')) returns contents of edit_th_winsize as a double

handles.th_winsize = get(handles.edit_th_winsize, 'String');
handles.th_winsize = str2num(handles.th_winsize);
startNdx = handles.startNdx;
endNdx = handles.endNdx;
handles = LoadFile(handles, handles.filenum);
handles.startNdx = startNdx;
handles.endNdx= endNdx;
handles = UpdateAxes(handles);

guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function edit_th_winsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_th_winsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit_th_winstep_Callback(hObject, eventdata, handles)
% hObject    handle to edit_th_winstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_th_winstep as text
%        str2double(get(hObject,'String')) returns contents of edit_th_winstep as a double
handles.th_winstep = get(handles.edit_th_winstep, 'String');
handles.th_winstep = str2num(handles.th_winstep);
startNdx = handles.startNdx;
endNdx = handles.endNdx;
handles = LoadFile(handles, handles.filenum);
handles.startNdx = startNdx;
handles.endNdx= endNdx;
handles = UpdateAxes(handles);

guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function edit_th_winstep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_th_winstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in check_UseRecWindow.
function check_UseRecWindow_Callback(hObject, eventdata, handles)
% hObject    handle to check_UseRecWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_UseRecWindow
end
