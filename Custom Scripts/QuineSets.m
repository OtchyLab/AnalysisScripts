function varargout = QuineSets(varargin)
% QUINESETS M-file for QuineSets.fig
%      QUINESETS, by itself, creates a new QUINESETS or raises the existing
%      singleton*.
%
%      H = QUINESETS returns the handle to a new QUINESETS or the handle to
%      the existing singleton*.
%
%      QUINESETS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QUINESETS.M with the given input arguments.
%
%      QUINESETS('Property','Value',...) creates a new QUINESETS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before QuineSets_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to QuineSets_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help QuineSets

% Last Modified by GUIDE v2.5 30-Jun-2010 19:41:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @QuineSets_OpeningFcn, ...
                   'gui_OutputFcn',  @QuineSets_OutputFcn, ...
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

% --- Executes just before QuineSets is made visible.
function QuineSets_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QuineSets (see VARARGIN)

% Choose default command line output for QuineSets
handles.output = hObject;

%Set startup salues for the GUI
set(handles.edit_MinPnts,'String',num2str(200));
set(handles.edit_Eps,'String',num2str(-1));
set(handles.edit_OrderedFileNum,'String',num2str(1));
set(handles.edit_Timescale,'String',num2str(0.1));
set(handles.edit_Chi,'String',num2str(0.01));
set(handles.edit_MaxClusSize,'String',num2str(25));

%Set startup values for the backend processing
handles.MinPnts = 200;
handles.Eps = -1;
handles.OrderedFileNum = 1;
handles.Timescale = 0.1;
handles.Chi = 0.01;
handles.MaxClusSize = 25;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes QuineSets wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = QuineSets_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

function [Eps] = FindEps(DataSet, MinPnts)
%Function attempts to estimate an optimal value for Eps from size of
%DataSet and neighborhood members threshold.

[m,n]=size(DataSet);

Eps=((prod(max(DataSet)-min(DataSet))*MinPnts*gamma(.5*n+1))/(m*sqrt(pi.^n))).^(1/n);

end %end function

function handles=UpdateAxesReachability(handles)
% Function updates GUI axesReachability display

%Plot OrderedFiles data on the Reachability plot
axes(handles.axesReachability);
hold off;
bar(handles.OrderedFiles(:,3), 'b');
xlim([1 length(handles.OrderedFiles)]);
ylim([0 round(max(handles.OrderedFiles(:,3))+1.5)]);
box on;
hold on;

%Plot a marker on the plot representing the selected syllable
line([handles.OrderedFileNum, handles.OrderedFileNum],[0,round(max(handles.OrderedFiles(:,3))+1.5)], 'color', 'r');
hold on;

%Plot shaded regions on the plot representing the edges of detected
%syllables
if isfield(handles,'ClusterSet');
    if ~isempty(handles.ClusterSet);
    syllStartPoints = handles.ClusterSet(:,1)';
    syllEndPoints = handles.ClusterSet(:,2)';
    if isempty(handles.ClusterSet)
        syllType = -1*ones(size(handles.ClusterSet,1),1);
    else
        syllType = -1*ones(size(handles.ClusterSet,1),1);
    end
    rects = [syllStartPoints; syllStartPoints; syllEndPoints; syllEndPoints; syllStartPoints];   
    color = zeros(3,length(syllStartPoints));
    color(3,syllType~=-1) = 1;
    color(1,syllType==-1) = 1;
%     if(handles.selectedSyll~=-1)
%         color(:,handles.selectedSyll) = [0,1,0];
%     end
    if ~isempty(syllStartPoints) && ~isempty(syllEndPoints)
        %draw nav bar rectangles.
        lims = ylim;
        y = [lims(1); lims(2); lims(2); lims(1); lims(1)];
        for(nSyll = 1:length(syllStartPoints))
            handles.rectHandles(nSyll,1) = patch(rects(:,nSyll),y, color(:,nSyll)', 'EdgeColor', 'none', 'FaceAlpha', .2);
            set(handles.rectHandles(nSyll,1), 'HitTest', 'off');
        end
        %set(handles.rectHandles(:,1), 'HitTest', 'off');
    end
    end
end

hold on;

set(handles.axesReachability, 'ButtonDownFcn', @cb_Reachability_click);

end

function handles=UpdateAxesFeatures(handles)
%Function updates GUI axesReachability display

%Plot Attributes data on the Features plot
axes(handles.axesFeatures);
hold off;
imagesc(handles.AttributeImage,[0 255]); colormap(jet);
xlim([1 length(handles.OrderedFiles)]);
ylim([.5 size(handles.AttributeImage,1)+.5]);
box on;
hold on;

%set(handles.axesFeatures, 'ButtonDownFcn', @cb_Features_click);

end

function cb_Reachability_click(hObject, evnt)
%Executes on clicking in axesReachability to set the current syllable

%Load handles so it's addressable in the function
handles = guidata(hObject);
mouseMode = get(get(hObject,'Parent'), 'SelectionType');
clickLocation = get(handles.axesReachability, 'CurrentPoint');
handles.OrderedFileNum = round(clickLocation(1,1));
set(handles.edit_OrderedFileNum, 'String',num2str(handles.OrderedFileNum));
handles.SampleSyl = handles.SubsetPointer(handles.OrderedFiles(handles.OrderedFileNum,1));
handles.AddMarker = 1;

handles = UpdateSylLabel(handles);
handles = UpdateSylIncExc(handles);

handles = UpdateAxesReachability(handles);
handles = UpdateAxesSampleSyl(handles);

guidata(hObject, handles);
end

function cb_templates_click(hObject, evnt)
%Executes on clicking in axesTemplate to remove syllable from the template

handles = guidata(hObject);
mouseMode = get(get(hObject,'Parent'), 'SelectionType');
clickLocation = get(handles.axesTemplate, 'CurrentPoint');
axes(handles.axesTemplate);

if(strcmp(mouseMode, 'open'))
    if(isfield(handles, 'template'))
        fs = 44100;
        edges = [0];
        allwavs = [];
        for(nWav = 1:length(handles.template.wavs))
            wav = handles.template.wavs(nWav).wav;
            label = handles.template.wavs(nWav).segType;
            allwavs = [allwavs; wav];
            edges = [edges, edges(end) + (length(wav) - 1)/fs];
        end
        nSelect = find(clickLocation(1,1) < edges);
        if(length(nSelect) > 0)
            nSelect = nSelect(1) - 1;
            %handles=updateTemplateFeatures(handles,nSelect,'delete');
            handles.template.wavs(nSelect) = [];
            handles = UpdateAxesTemplate(handles);
        end
    end
end

guidata(hObject, handles);
end

function handles=UpdateAxesSampleSyl(handles)
%Function updates GUI axesReachability display

%Plot OrderedFiles data on the Reachability plot
axes(handles.axesSampleSyl);
hold off;

%Location set by GUI
if isfield(handles,'WavFileLoc')
    SylFolder = handles.WavFileLoc;
else
    SylFolder = 'C:\Users\Tim\Desktop\Test Data for Quine\Syllables Exported';
end

%HardCoded for now
SylLoc = [SylFolder '\' char(handles.FeatureKeys(handles.SampleSyl))];
[handles.wav handles.fs] = wavread(SylLoc);
displaySpecgramQuick(handles.wav, handles.fs, [0,10000]);
xlabel('');
ylabel('');
xlim([0 handles.Timescale]);
set(handles.axesSampleSyl,'XTick',[]);
hold on;

end

function [handles] = UpdateAxesTemplate(handles)
%Updates the GUI axesTemplate display

%Clear the current axes
cla(handles.axesTemplate);

%Only execute if a template exists/is loaded
if(isfield(handles, 'template'))
    fs = 44100;
    allwavs = [];
    edges = [0];
    for(nWav = 1:length(handles.template.wavs))
        wav = handles.template.wavs(nWav).wav;
        label = handles.template.wavs(nWav).segType;
        allwavs = [allwavs; wav];
        edges = [edges, edges(end) + (length(wav) - 1)/fs];
    end
    
    %Add additional buffer to the end of the display file
    if(edges(end) < 1)
        allwavs(end+1:fs) = 0;
    end
    
    %Plot Spectro on the axes
    axes(handles.axesTemplate);
    displaySpecgramQuick(allwavs, fs, [0,10000]);
    
    %Individuate the syls and label with the Syllable type
    for(nWav = 1:length(edges)-1)
        l = line([edges(nWav),edges(nWav)], ylim);
        set(l,'Color','k');
        set(l,'LineWidth', 2);
        mid = mean(edges([nWav,nWav+1]));
        if(handles.template.wavs(nWav).segType ~= -1)
            t = text(mid, mean(ylim), num2str(handles.template.wavs(nWav).segType), 'FontWeight', 'bold');
        end
    end
    set(handles.axesTemplate, 'ButtonDownFcn', @cb_templates_click);
end

end

function [handles] = UpdateOrdering(handles, DataSet)
%Handles the sequence for OPTICS sorting the syllables and updating the GUI
%displays accordingly.  The output replaces the current
%handles.OrderedFiles structure.

%PCA the passed DataSet
[PCA_set] = getPCA(DataSet);
if handles.Eps == -1
    handles.Eps = FindEps(PCA_set,handles.MinPnts);
    set(handles.text_EpsCalc,'String',num2str(handles.Eps));
end

%Order Dataset
[handles.OrderedFiles, handles.Eps]=OPTICS(PCA_set, handles.Eps, handles.MinPnts);

%Temp for Debugging:
%load('Backup Post Ordering.mat')
%handles.OrderedFiles=OrderedFiles;

%Plot the OrderedFiles data to the GUI
handles = UpdateAxesReachability(handles);
%Create discretized attribute image for the entire dataset (hence, the [])
[handles.AttributeImage] = ExtractAttributePlots(handles.OrderedFiles, DataSet, []);
%Plot the OrderedFiles data to the GUI
handles = UpdateAxesFeatures(handles);
%Set flag to allow clustering
handles.IsOrdered = 1;
set(handles.text_NumFiles,'String',num2str(length(handles.OrderedFiles)));

end

function [handles, ClustMem] = WhichCluster(handles, Target)
%Determines which found cluster the selected syllable belongs to.  If
%multiple clusters found, it asks the user which to select as the
%appropriate one.

Candidates = [];
for i=1:size(handles.ClusterSet,1)
    Clust = handles.SubsetPointer(handles.OrderedFiles(handles.ClusterSet(i,1):handles.ClusterSet(i,2),1));
    if ismember(Target,Clust)
        Candidates = [Candidates; i];
    end
end
if isempty(Candidates)
    error('Not a member of any defined cluster.')
    return
elseif length(Candidates)==2
    BT1 = [num2str(handles.ClusterSet(Candidates(1),1)) ' to ' num2str(handles.ClusterSet(Candidates(1),2))];
    BT2 = [num2str(handles.ClusterSet(Candidates(2),1)) ' to ' num2str(handles.ClusterSet(Candidates(2),2))];
    BT4 = 'Cancel';
    WhichClust = questdlg('Syllable found to be a part of more than one identified cluster.  Which to used?', 'WARNING!', BT1, BT2, BT4);
elseif length(Candidates)==3
    BT1 = [num2str(handles.ClusterSet(Candidates(1),1)) ' to ' num2str(handles.ClusterSet(Candidates(1),2))];
    BT2 = [num2str(handles.ClusterSet(Candidates(2),1)) ' to ' num2str(handles.ClusterSet(Candidates(2),2))];
    BT3 = [num2str(handles.ClusterSet(Candidates(3),1)) ' to ' num2str(handles.ClusterSet(Candidates(3),2))];
    BT4 = 'Cancel';
    WhichClust = questdlg('Syllable found to be a part of more than one identified cluster.  Which to used?', 'WARNING!', BT1, BT2, BT3, BT4);
elseif length(Candidates)==1
    BT1 = [num2str(handles.ClusterSet(Candidates(1),1)) ' to ' num2str(handles.ClusterSet(Candidates(1),2))];
    WhichClust = BT1;
end

switch WhichClust
    case BT1
        ClustMem = handles.SubsetPointer(handles.OrderedFiles(handles.ClusterSet(Candidates(1),1):handles.ClusterSet(Candidates(1),2),1));
    case BT2
        ClustMem = handles.SubsetPointer(handles.OrderedFiles(handles.ClusterSet(Candidates(2),1):handles.ClusterSet(Candidates(2),2),1));
    case BT3
        ClustMem = handles.SubsetPointer(handles.OrderedFiles(handles.ClusterSet(Candidates(3),1):handles.ClusterSet(Candidates(3),2),1));
    case BT4
        return
end
end

 function handles = UpdateSylLabel(handles)
%This function updates the syllable label box on the GUI

if handles.SylTypes(handles.SampleSyl) ~= -1
    set(handles.edit_SylLabel,'String',num2str(handles.SylTypes(handles.SampleSyl)));
else
    set(handles.edit_SylLabel,'String','');
end

 end

 function handles = UpdateSylIncExc(handles)
%This function updates the inclusion/exclusion text on the GUI
 
if handles.SubSet(handles.SampleSyl)== -1
    set(handles.text_SylIncExc, 'ForegroundColor','r'); 
    set(handles.text_SylIncExc, 'String','Excluded');
else
    set(handles.text_SylIncExc, 'ForegroundColor','g'); 
    set(handles.text_SylIncExc, 'String','Included');
end

 end
 
 function [SimAnal] = SimAnalysis(handles)
 %Function take the Features and Labeled Syls structures (in their current
 %state) and calculated similarity, variability, etc. for each ID'd
 %syllable type.  The SubSet structure is completely ignored.

SimAnal =[];
FeatureSet = handles.SimFeatures; %Can be changed to handles.FeatureFile to process 33-param file
SylLabels = handles.SylTypes;

SimAnal.UniqueSyls = sort(unique(SylLabels));

%Remove syllables to be ignored
SimAnal.UniqueSyls = SimAnal.UniqueSyls(SimAnal.UniqueSyls>0 & SimAnal.UniqueSyls<100);
SimAnal.SimSet = FeatureSet(SylLabels>0 & SylLabels<100,:);
SimAnal.SylSet = SylLabels(SylLabels>0 & SylLabels<100);

%Find number of syls to calc stats for
NumSyls = length(SimAnal.UniqueSyls); 

%Z-score the matrix of ID'd features
[r c] = size(SimAnal.SimSet);
m_SimSet = mean(SimAnal.SimSet);
std_SimSet = std(SimAnal.SimSet);
for i=1:c
    Zscore(:,i) = (SimAnal.SimSet(:,i)-m_SimSet(i))./std_SimSet(i);
end
SimAnal.SimSet = Zscore;

%Calculates the maximum difference (i.e., distance) between two ID'd syls
SimAnal.maxDist = 3.5; %max(pdist(SimAnal.SimSet));  %running about 10.8


SimStack = [];
TotalSyls = 0;

for i=1:NumSyls
    SimAnal.Stats(i).Dist = pdist(SimAnal.SimSet(SimAnal.SylSet==SimAnal.UniqueSyls(i)));
    SimAnal.Stats(i).Sim = (SimAnal.maxDist-SimAnal.Stats(i).Dist).*100./SimAnal.maxDist;
    SimAnal.Stats(i).meanSim = mean(SimAnal.Stats(i).Sim);
    SimAnal.Stats(i).Var = (100-SimAnal.Stats(i).meanSim)*100/50;
    SimAnal.Stats(i).Centroid = mean(SimAnal.SimSet(SimAnal.SylSet==SimAnal.UniqueSyls(i),:));
    SimAnal.Stats(i).meanDist = mean(SimAnal.Stats(i).Dist);
    SimAnal.Stats(i).stdDist = std(SimAnal.Stats(i).Dist);
    SimAnal.Stats(i).numSyl = length(SimAnal.SimSet(SimAnal.SylSet==SimAnal.UniqueSyls(i)));
    [m,n]=size(SimAnal.SimSet(SimAnal.SylSet==SimAnal.UniqueSyls(i),:));
    SimAnal.Stats(i).Radii=sqrt(sum((((ones(m,1)*SimAnal.Stats(i).Centroid)-SimAnal.SimSet(SimAnal.SylSet==SimAnal.UniqueSyls(i),:)).^2)'));
    if n==1
        SimAnal.Stats(i).Radii=abs((ones(m,1)*i-SimAnal.SimSet(SimAnal.SylSet==SimAnal.UniqueSyls(i),:)))';
    end
    SimStack = [SimStack;SimAnal.Stats(i).Sim'];
    TotalSyls = TotalSyls + SimAnal.Stats(i).numSyl;
end

SimAnal.SetVar = (100-mean(SimStack))*100/50;

%Display the result on a new figure

%Determine number of figures required
if (NumSyls/4)>round(NumSyls/4)
    NumFigs = round(NumSyls/4)+1;
else
    NumFigs = round(NumSyls/4);
    MkNew = 1;
end

CurrSyl = 1;
max_aax=0;
max_aay=0;
max_bbx=0;
max_bby=0;

for fig = 1:NumFigs
    figure;
    i=1;
    while (i < 5) && (CurrSyl < NumSyls+1)
        aa(CurrSyl)=subplot(4,3,(3*i-2));
        hist(SimAnal.Stats(CurrSyl).Sim,50);
        max_aax = max(xlim,max_aax);
        max_aay = max(ylim,max_aay);
        title(['Similarity distribution for Syllable ' num2str(SimAnal.UniqueSyls(CurrSyl))])
        box on
        %hold on
        bb(CurrSyl)=subplot(4,3,(3*i-1));
        hist(SimAnal.Stats(CurrSyl).Radii,50);
        max_bbx = max(xlim,max_bbx);
        max_bby = max(ylim,max_bby);
        %title(['Radii distribution for Syllable ' num2str(SimAnal.UniqueSyls(i))])
        box on
        %hold on
        ah=subplot(4,3,(3*i));
        disp_string=['Mean Similarity: ' num2str(SimAnal.Stats(CurrSyl).meanSim) sprintf('\n') 'Variability: ' num2str(SimAnal.Stats(CurrSyl).Var) sprintf('\n') '# of Syllables: ' num2str(SimAnal.Stats(CurrSyl).numSyl)];
        text(0, 0.25, disp_string);
        set(ah,'Visible','off');
        %title(['Stats for Syllable ' num2str(SimAnal.UniqueSyls(i))])
        %box on
        %hold on
        i=i+1;
        CurrSyl = CurrSyl+1;
    end
end

set(aa(:),'xlim',[50, max_aax(2)]);
set(bb(:),'xlim',[0 10]);
set(aa(:),'ylim',max_aay);
set(bb(:),'ylim',max_bby);

if exist('MkNew','var')
    figure;
    i=2;
end
ah=subplot(4,3,(3*(i)-2:3*(i)));
disp_string=['Set Variability: ' num2str(SimAnal.SetVar) sprintf('\n') 'Total Number of Syllables Analyzed: ' num2str(TotalSyls)];
text(.25,.25,disp_string);
set(ah,'Visible','off');
%title('Stats for Syllable Set')
box on
hold on

 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The below control the GUI buttons, boxes, etc.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in push_PrevFile.
function push_PrevFile_Callback(hObject, eventdata, handles)
% hObject    handle to push_PrevFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.OrderedFileNum < 2
    error('Requested file number out of range.')
    return
else
    handles.OrderedFileNum = handles.OrderedFileNum - 1;
    set(handles.edit_OrderedFileNum, 'String', num2str(handles.OrderedFileNum))
    handles.SampleSyl = handles.SubsetPointer(handles.OrderedFiles(handles.OrderedFileNum,1));
    handles = UpdateSylLabel(handles);
    handles = UpdateSylIncExc(handles);
    handles = UpdateAxesReachability(handles);
    handles = UpdateAxesSampleSyl(handles);
end

guidata(hObject, handles);
end

% --- Executes on button press in push_NextFile.
function push_NextFile_Callback(hObject, eventdata, handles)
% hObject    handle to push_NextFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.OrderedFileNum > (length(handles.OrderedFiles)-1)
    error('Requested file number out of range.')
    return
else
    handles.OrderedFileNum = handles.OrderedFileNum + 1;
    set(handles.edit_OrderedFileNum, 'String', num2str(handles.OrderedFileNum))
    handles.SampleSyl = handles.SubsetPointer(handles.OrderedFiles(handles.OrderedFileNum,1));
    handles = UpdateSylLabel(handles);
    handles = UpdateSylIncExc(handles);
    handles = UpdateAxesReachability(handles);
    handles = UpdateAxesSampleSyl(handles);
end

guidata(hObject, handles);

end

function edit_OrderedFileNum_Callback(hObject, eventdata, handles)
% hObject    handle to edit_OrderedFileNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_OrderedFileNum as text
%        str2double(get(hObject,'String')) returns contents of edit_OrderedFileNum as a double
handles.OrderedFileNum = get(handles.edit_OrderedFileNum, 'String');
handles.OrderedFileNum = str2num(handles.OrderedFileNum);

if handles.OrderedFileNum == 0 || handles.OrderedFileNum>length(handles.OrderedFiles)
    error('Requested file number out of range.')
    return
else
    handles.SampleSyl = handles.SubsetPointer(handles.OrderedFiles(handles.OrderedFileNum,1));
    handles = UpdateSylLabel(handles);
    handles = UpdateSylIncExc(handles);
    
    handles = UpdateAxesReachability(handles);
    handles = UpdateAxesSampleSyl(handles);
end

guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function edit_OrderedFileNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_OrderedFileNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in check_NameAllInCluster.
function check_NameAllInCluster_Callback(hObject, eventdata, handles)
% hObject    handle to check_NameAllInCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_NameAllInCluster
end

function edit_SylLabel_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SylLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SylLabel as text
%        str2double(get(hObject,'String')) returns contents of edit_SylLabel as a double
handles.SylLabel = get(handles.edit_SylLabel, 'String');
handles.SylLabel = str2num(handles.SylLabel);

guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function edit_SylLabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SylLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in push_ApplySylLabel.
function push_ApplySylLabel_Callback(hObject, eventdata, handles)
% hObject    handle to push_ApplySylLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.check_NameAllInCluster,'Value') == 0
    handles.SylTypes(handles.SampleSyl) = handles.SylLabel;
else
    [handles, ClustMem] = WhichCluster(handles, handles.SampleSyl);
    handles.SylTypes(ClustMem) = handles.SylLabel;
end

guidata(hObject, handles);

end

% --- Executes on button press in push_OpenFeatureFile.
function push_OpenFeatureFile_Callback(hObject, eventdata, handles)
% hObject    handle to push_OpenFeatureFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Select the .mat file containing the Feature Set
%handles.FeatFileLoc = uigetfile('*.mat', 'Features File');
[handles.FeatureFileName, handles.FeatureFileLoc]=uigetfile('*.mat', 'Features File');
load([handles.FeatureFileLoc '\' handles.FeatureFileName]);

%This line will process the 33-param file
handles.FeatureFile = FeatureSet;

%This line will process the 10-param file
handles.SimFeatures = FeatureTotals;

%Contains the file names of each syllable that corresponds to a row in
%handles.FeatureSet
handles.FeatureKeys = FeatureKeys;

%Copy and/or initialize Syllable Type and SubSet arrays
if ~exist('SylTypes')
    handles.SylTypes = -1*ones(length(handles.FeatureFile),1);
else
    handles.SylTypes = SylTypes;
end

handles.SubSet = (1:length(handles.FeatureFile))';

%Reset ordering and clustering flags
handles.IsOrdered = 0;
handles.IsClustered = 0;

guidata(hObject, handles);

end

% --- Executes on button press in push_OrderData.
function push_OrderData_Callback(hObject, eventdata, handles)
% hObject    handle to push_OrderData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Check to see that all the necessary info is loaded
if isempty(handles.MinPnts) || handles.MinPnts<=1
    error('Enter minimum points to define neightborhood greater than 1')
    return
elseif isempty(handles.Eps) || (handles.Eps<=0 && handles.Eps~=-1)
    error('Enter minimum points to define neightborhood greater than 0 or enter -1 to calculate')
    return
elseif isempty(handles.FeatureFile)
    error('Must load Feature File before sorting')
    return
else
    %Confirm that you want to overwrite OrderedFiles, IF it already exists
    if isfield(handles, 'OrderedFiles')
        Confirm = questdlg('Continuing will overwrite current Ordered List.  Continue?', 'WARNING!', 'No Problem', 'Oh, Shit!', 'Cancel');
    else
        Confirm = 'No Problem';
    end
    switch Confirm
        case 'No Problem'
            %Create the 
            handles.SubsetPointer = handles.SubSet(handles.SubSet~=-1);
            handles.SubsetFeatureFile = handles.FeatureFile(handles.SubsetPointer,:);
            
            handles.ClusterSet = [];
            
            [r c] = size(handles.SubsetFeatureFile);
            
            %Z-score the entire matrix
            m_FeatureSet = mean(handles.SubsetFeatureFile);
            std_FeatureSet = std(handles.SubsetFeatureFile);
            for i=1:c
                Zscore(:,i) = (handles.SubsetFeatureFile(:,i)-m_FeatureSet(i))./std_FeatureSet(i);
            end
            handles.z_SubsetFeatures = Zscore;
            [handles] = UpdateOrdering(handles, handles.z_SubsetFeatures);
            
        case 'Oh, Shit!'
           
        case 'Cancel'
            
    end
end
guidata(hObject, handles);

end

function edit_MinPnts_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MinPnts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MinPnts as text
%        str2double(get(hObject,'String')) returns contents of edit_MinPnts as a double
handles.MinPnts = get(handles.edit_MinPnts, 'String');
handles.MinPnts = str2num(handles.MinPnts);

guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function edit_MinPnts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MinPnts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit_Eps_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Eps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Eps as text
%        str2double(get(hObject,'String')) returns contents of edit_Eps as a double
handles.Eps = get(handles.edit_Eps, 'String');
handles.Eps = str2num(handles.Eps);

guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function edit_Eps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Eps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit_Chi_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Chi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Chi as text
%        str2double(get(hObject,'String')) returns contents of edit_Chi as a double
handles.Chi = get(handles.edit_Chi, 'String');
handles.Chi = str2num(handles.Chi);

guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function edit_Chi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Chi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit_MaxClusSize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MaxClusSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MaxClusSize as text
%        str2double(get(hObject,'String')) returns contents of edit_MaxClusSize as a double
handles.MaxClusSize = get(handles.edit_MaxClusSize, 'String');
handles.MaxClusSize = str2num(handles.MaxClusSize);

guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function edit_MaxClusSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MaxClusSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in push_ClusterData.
function push_ClusterData_Callback(hObject, eventdata, handles)
% hObject    handle to push_ClusterData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.IsOrdered ~= 1;%Check that new data has been loaded and ordered for clustering
    error('You must load and order a dataset before clustering')
    return
elseif isempty(handles.MaxClusSize) || (handles.MaxClusSize<=0 || handles.MaxClusSize>100)
    error('Enter Max Cluster Size value greater than 0 and less than 100')
    return
elseif isempty(handles.Chi) || (handles.Chi<=0 || handles.Chi>100)
    error('Enter Chi-steepness value greater than 0 and less than 100')
    return
else
    %Identify Chi-steep regions from OrderedFiles
    [SDASet, SUASet, IndexIsDownUp] = SteepRegions(handles.OrderedFiles, handles.Chi, handles.MinPnts);

    %Extract clusters fromt the set of Chi-steep regions
    [handles.ClusterSetRaw] = ExtractOPTICSClusters(handles.OrderedFiles, SDASet, SUASet, handles.Chi, handles.MinPnts);

    %Filter the found set of clusters by size --> only keep those with smaller
    %than the max % cluster size set in the GUI
    ClusDiff = handles.ClusterSetRaw(:,2)-handles.ClusterSetRaw(:,1);
    handles.ClusterSet = handles.ClusterSetRaw(ClusDiff<((handles.MaxClusSize/100)*length(handles.OrderedFiles)),:);
    set(handles.text_NumClust,'String',num2str(length(handles.ClusterSet)));
    
    %Set the clustered flag to allow future functions
    if isempty(handles.ClusterSet)
        handles.IsClustered = 0;
    else
        handles.IsClustered = 1;
        %Update Reachability plot to overlay the found clusters
        handles = UpdateAxesReachability(handles);
    end
end

guidata(hObject, handles);

end

% --- Executes on button press in push_ReloadOrderedSet.
function push_ReloadOrderedSet_Callback(hObject, eventdata, handles)
% hObject    handle to push_ReloadOrderedSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Temp for Debugging:
%load('OrderedFiles.mat')

%Select the .mat file containing the Ordered Set
[handles.OrderFileName, handles.OrderFileLoc]=uigetfile('*.mat', 'Ordering File');
load([handles.OrderFileLoc '\' handles.OrderFileName]);

%Load Fields to variables
handles.FeatureFile = FeatureFile; handles.FeatureKeys = FeatureKeys; handles.SylTypes=SylTypes;
handles.SubSet = SubSet; handles.OrderedFiles = OrderedFiles; handles.IsOrdered = 1; handles.IsClustered = 0;
handles.SubsetPointer = SubsetPointer; handles.SubsetFeatureFile = SubsetFeatureFile; 

if exist('SimFeatures') == 1
    handles.SimFeatures = SimFeatures;
end

%Plot the OrderedFiles data to the GUI
handles = UpdateAxesReachability(handles);
%Create discretized attribute image for the entire dataset (hence, the [])
[handles.AttributeImage] = ExtractAttributePlots(handles.OrderedFiles, handles.SubsetFeatureFile, []);
%Plot the OrderedFiles data to the GUI
handles = UpdateAxesFeatures(handles);
set(handles.text_NumFiles,'String',num2str(length(handles.OrderedFiles)));

guidata(hObject, handles);

end

function edit_Timescale_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Timescale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Timescale as text
%        str2double(get(hObject,'String')) returns contents of edit_Timescale as a double

handles.Timescale = get(handles.edit_Timescale, 'String');
handles.Timescale = str2num(handles.Timescale);

if ~isempty(handles.SampleSyl)
    handles = UpdateAxesSampleSyl(handles);
end

guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function edit_Timescale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Timescale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in push_PlotProjection.
function push_PlotProjection_Callback(hObject, eventdata, handles)
% hObject    handle to push_PlotProjection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure;
if get(handles.popup_Proj1, 'Value') == 1
    ax1 = 1:length(handles.OrderedFiles);
else
    FeatSel = get(handles.popup_Proj1, 'Value')-1;
    ax1 = handles.FeatureFile(handles.SubsetPointer,FeatSel);
end

if get(handles.popup_Proj2, 'Value') == 1
    ax2 = 1:length(handles.OrderedFiles);
else
    FeatSel = get(handles.popup_Proj2, 'Value')-1;
    ax2 = handles.FeatureFile(handles.SubsetPointer,FeatSel);
end
scatter(ax1,ax2,'.');
xlim([min(ax1)-(0.1*(max(ax1)-min(ax1))), max(ax1)+(0.1*(max(ax1)-min(ax1)))]);
ylim([min(ax2)-(0.1*(max(ax2)-min(ax2))), max(ax2)+(0.1*(max(ax2)-min(ax2)))]);

end

% --- Executes on button press in push_SaveOrderedSet.
function push_SaveOrderedSet_Callback(hObject, eventdata, handles)
% hObject    handle to push_SaveOrderedSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Request user input for where to save Ordered Set file
[handles.OrderedSetFileName, handles.OrderedSetFileLoc]=uiputfile('features.mat','Save OrderedSet as:');

%Load Fields to variables
FeatureFile = handles.FeatureFile; FeatureKeys = handles.FeatureKeys; SylTypes=handles.SylTypes;
SubSet = handles.SubSet; OrderedFiles = handles.OrderedFiles; IsOrdered = handles.IsOrdered;
SubsetPointer = handles.SubsetPointer; SubsetFeatureFile=handles.SubsetFeatureFile; SimFeatures = handles.SimFeatures;
            
%Save the files necessary to restore the session
save([handles.OrderedSetFileLoc '\' handles.OrderedSetFileName], 'FeatureFile', 'FeatureKeys', 'SylTypes', 'SubSet','OrderedFiles','IsOrdered','SubsetPointer', 'SubsetFeatureFile','SimFeatures');

guidata(hObject, handles);

end

% --- Executes on selection change in popup_Proj1.
function popup_Proj1_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Proj1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Proj1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Proj1
end

% --- Executes during object creation, after setting all properties.
function popup_Proj1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Proj1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in popup_Proj2.
function popup_Proj2_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Proj2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Proj2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Proj2
end

% --- Executes during object creation, after setting all properties.
function popup_Proj2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Proj2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in check_MarkClust.
function check_MarkClust_Callback(hObject, eventdata, handles)
% hObject    handle to check_MarkClust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_MarkClust
end


% --- Executes on button press in push_OpenTemplate.
function push_OpenTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to push_OpenTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'template')
    answer = questdlg('This will overwrite the current template.  Continue?', 'WARNING!', 'Yes', 'No', 'No');
else
    answer = 'Yes';
end

switch answer
    case 'Yes'
        uiload;
        handles.template = templates;
        handles = UpdateAxesTemplate(handles);
    case 'No' 
end

guidata(hObject, handles);

end

% --- Executes on button press in push_SaveTemplate.
function push_SaveTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to push_SaveTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'template')
    %Request user input for where to save Ordered Set file
    [TemplateFileName, TemplateFileLoc]=uiputfile('Template.mat','Save Template as:');

    %Load Fields to variables
    templates = handles.template;
    templateFeatures =[]; %This is empty and required (?) for compatability with legacy programs
            
    %Save the files necessary to restore the session
    save([TemplateFileLoc '\' TemplateFileName], 'templates', 'templateFeatures');
else
    warndlg('No templates to save.');
    uiwait;
end
end

% --- Executes on button press in push_AddTemplate.
function push_AddTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to push_AddTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.SylTypes(handles.SampleSyl) ~=-1
    if isfield(handles, 'template')
        NewSylSpot = length(handles.template.wavs)+1;
        handles.template.wavs(NewSylSpot).wav = handles.wav;
        handles.template.wavs(NewSylSpot).segType = handles.SylTypes(handles.SampleSyl);
        handles.template.wavs(NewSylSpot).features = handles.FeatureFile(handles.SampleSyl);
        handles.template.wavs(NewSylSpot).fs = 44100;
        handles = UpdateAxesTemplate(handles);
    else
        handles.template.wavs.wav = [];
        handles.template.wavs(1).wav = handles.wav;
        handles.template.wavs.segType = [];
        handles.template.wavs(1).segType = handles.SylTypes(handles.SampleSyl);
        handles.template.wavs.features = [];
        handles.template.wavs(1).features = handles.FeatureFile(handles.SampleSyl,:);
        handles.template.wavs.fs = [];
        handles.template.wavs(1).fs = 44100;
        handles = UpdateAxesTemplate(handles);
    end
else
    error('You must label the syllable before adding to the template')
    return
end

guidata(hObject, handles);

end

% --- Executes on button press in check_NoSylFiles.
function check_NoSylFiles_Callback(hObject, eventdata, handles)
% hObject    handle to check_NoSylFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_NoSylFiles
end

% --- Executes on button press in push_SetWavLoc.
function push_SetWavLoc_Callback(hObject, eventdata, handles)
% hObject    handle to push_SetWavLoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Request user input for where to find wav-files for SylableSample Preview
handles.WavFileLoc = uigetdir(pwd,'Location of Syllable or Song Wav Files');

guidata(hObject, handles);

end

% --- Executes on button press in push_RemoveSyl.
function push_RemoveSyl_Callback(hObject, eventdata, handles)
% hObject    handle to push_RemoveSyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.SubSet(handles.SampleSyl) = -1;
handles = UpdateSylIncExc(handles);

guidata(hObject, handles);

end

% --- Executes on button press in push_RemoveClust.
function push_RemoveClust_Callback(hObject, eventdata, handles)
% hObject    handle to push_RemoveClust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles, Clust] = WhichCluster(handles, handles.SampleSyl);

handles.SubSet(Clust) = -1;
handles = UpdateSylIncExc(handles);

guidata(hObject, handles);

end

% --- Executes on button press in push_AddSyl.
function push_AddSyl_Callback(hObject, eventdata, handles)
% hObject    handle to push_AddSyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.SubSet(handles.SampleSyl) = handles.SampleSyl;
handles = UpdateSylIncExc(handles);

guidata(hObject, handles);

end


% --- Executes on button press in push_AddCluster.
function push_AddCluster_Callback(hObject, eventdata, handles)
% hObject    handle to push_AddCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles, Clust] = WhichCluster(handles, handles.SampleSyl);

handles.SubSet(Clust) = Clust;
handles = UpdateSylIncExc(handles);

guidata(hObject, handles);
end


% --- Executes on button press in push_AddAll.
function push_AddAll_Callback(hObject, eventdata, handles)
% hObject    handle to push_AddAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Confirm = questdlg('Continuing will reset all selections since last Ordering.  Continue?', 'WARNING!', 'No Problem', 'Oh, Shit!', 'Cancel');

switch Confirm
    case 'No Problem'
        handles.SubSet = 1:length(handles.OrderedFiles);
        handles = UpdateSylIncExc(handles);
        handles = UpdateAxesSampleSyl(handles);
    otherwise
end

guidata(hObject, handles);
end


% --- Executes on button press in push_RemoveAllLabeled.
function push_RemoveAllLabeled_Callback(hObject, eventdata, handles)
% hObject    handle to push_RemoveAllLabeled (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Confirm = questdlg('Continuing will overwrite all selections since last Ordering.  Continue?', 'WARNING!', 'No Problem', 'Oh, Shit!', 'Cancel');

switch Confirm
    case 'No Problem'
        %Find all labelled syls in current set
        LabeledSyls = handles.SubsetPointer(handles.SylTypes(handles.SubsetPointer)~=-1);
        handles.SubSet(LabeledSyls) = -1;
        
        handles = UpdateSylIncExc(handles);
        handles = UpdateAxesSampleSyl(handles);
    otherwise
end

guidata(hObject, handles);
end

% --- Executes on button press in push_SimAnal.
function push_SimAnal_Callback(hObject, eventdata, handles)
% hObject    handle to push_SimAnal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.SimAnal] = SimAnalysis(handles);

guidata(hObject, handles);
end

% --- Executes on button press in push_ExportSimAnal.
function push_ExportSimAnal_Callback(hObject, eventdata, handles)
% hObject    handle to push_ExportSimAnal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles, 'SimAnal')
    %Request user input for where to save Ordered Set file
    [SimFileName, SimFileLoc]=uiputfile('Similarity.mat','Save Similarity Data as:');
    
    %Add template
    if isfield(handles, 'SimAnal.Template')
        handles.SimAnal.Template = handles.template;
    end
    
    %Load Fields to variables
    SimAnal = handles.SimAnal;
    FeatureSet = handles.FeatureFile;
    FeatureFile = handles.SimFeatures;
    SylLables = handles.SylTypes;
    FeatureKeys = handles.FeatureKeys;
            
    %Save the files necessary to restore the session
    save([SimFileLoc '\' SimFileName], 'SimAnal', 'FeatureSet', 'FeatureFile', 'SylLables', 'FeatureKeys');
else
    warndlg('No Analysis files to save.');
    uiwait;
end

end

% --- Executes on button press in push_OpenAnnotation.
function push_OpenAnnotation_Callback(hObject, eventdata, handles)
% hObject    handle to push_OpenAnnotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on button press in push_SaveAnnotation.
function push_SaveAnnotation_Callback(hObject, eventdata, handles)
% hObject    handle to push_SaveAnnotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end
