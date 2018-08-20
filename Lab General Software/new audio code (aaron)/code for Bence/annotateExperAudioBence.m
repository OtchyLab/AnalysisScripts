function varargout = annotateExperAudioBence(varargin)
% ANNOTATEEXPERAUDIOBENCE M-file for annotateExperAudioBence.fig
%      ANNOTATEEXPERAUDIOBENCE, by itself, creates a new ANNOTATEEXPERAUDIOBENCE or raises the existing
%      singleton*.
%
%      H = ANNOTATEEXPERAUDIOBENCE returns the handle to a new ANNOTATEEXPERAUDIOBENCE or the handle to
%      the existing singleton*.
%
%      ANNOTATEEXPERAUDIOBENCE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANNOTATEEXPERAUDIOBENCE.M with the given input arguments.
%
%      ANNOTATEEXPERAUDIOBENCE('Property','Value',...) creates a new ANNOTATEEXPERAUDIOBENCE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before annotateExperAudioBence_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to annotateExperAudioBence_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help annotateExperAudioBence


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @annotateExperAudioBence_OpeningFcn, ...
                   'gui_OutputFcn',  @annotateExperAudioBence_OutputFcn, ...
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


% --- Executes just before annotateExperAudioBence is made visible.
function annotateExperAudioBence_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to annotateExperAudioBence (see VARARGIN)

% Choose default command line output for annotateExperAudioBence
handles.output = hObject;
handles.exper = [];
handles.filenum = -1;
handles.selectedSyll = -1;
handles.startNdx = -1;
handles.endNdx = -1;
handles.txtHandles = [];
handles.rectHandles = [];
handles.navRect = [];
handles.spike_filename=[];
handles.cell_no=1; 
handles.audioAnnotation = mhashtable;
handles.bWaitingForAddClick = false;
handles.bWaitingForSeparationClick = false;
handles.threshold=0.8;
handles.edges=0.6;
handles.drugindex=1;
handles.spikemarker=[];
set(handles.popupDrug,'Value',handles.drugindex);
contents = get(handles.popupDrug,'String');
handles.drugstatus=contents{handles.drugindex};
set(handles.editPowerThreshold,'String',num2str(handles.threshold));
set(handles.editEdgeThreshold,'String',num2str(handles.edges));
set(handles.checkboxShowSpikes,'Value',0);
set(handles.checkboxShowSpikes,'Value',0);
set(handles.figure1,'interruptible','off');
set(handles.axesNavBar,'interruptible','off');
% 
% buttonLoadExper
% editFileNum
% buttonPrevFile
% buttonNextFile
% axesSpecgram

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes annotateExperAudioBence wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = annotateExperAudioBence_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if(~isequal(handles.exper,[]))
    saveExper(handles.exper);
end
varargout{1} = handles.output;


% --- Executes on button press in buttonNewAnnotation.
function buttonNewAnnotation_Callback(hObject, eventdata, handles)
% hObject    handle to buttonNewAnnotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%If annotation is already open check whether they want to append or clear.
if(isfield(handles, 'annotationFileName'))
    ansButton = questdlg('Would you like to clear the current annotations, or append to them?','',{'Clear','Append'});
end

%Get the new filename.
[file,path] = uiputfile('*.mat', 'Create a .mat file for the audioAnnotation:');

if isequal(file,0) || isequal(path,0)
else
    %if they don't his cancel..
    if(isfield(handles,'annotationFileName') && strcmp(ansButton, 'Clear'))
        %already have an annotation object loaded, then save clear current annotation.
        saveAnnotation(handles);
        guidata(hObject, handles);
        buttonClearAll_Callback(hObject, eventdata, handles);
        handles = guidata(hObject);
        handles.audioAnnotation.delete;        
    elseif(~isfield(handles,'annotationFilename'))
        handles.audioAnnotation = mhashtable;
    end
    handles.annotationFileName = [path,filesep,file];
    saveAnnotation(handles);
end
guidata(hObject, handles);

function saveAnnotation(handles)
if(isfield(handles,'annotationFileName'))
    aaSaveHashtable(handles.annotationFileName, handles.audioAnnotation);
else
    warndlg('You have not yet created an annotation.');
    uiwait;
end

function handles=saveCell(handles)

if(isfield(handles,'cell_allrecs'))
    presentCell=handles.cell_allrecs; 
     %if there are spikes for that record
    if (~isfield(handles,'sorted_spikes') | isempty(handles.sorted_spikes))
       % if there are no spikes for this record don't save anything
    else
        %insert the present record and then sort the structure based on filenumber before saving   
        handles.cell.spikes=handles.sorted_spikes;
        handles.cell.noise=handles.unsorted_spikes;
        filenumVector=getFieldVector(presentCell,'filenum'); %get the filenums in the present cell
        index=find(filenumVector==handles.filenum);
        if (isempty(index)) %if there is no record for the present file number then insert it appropriately
            index=find(filenumVector>handles.filenum,1,'first');
            if (isempty(index))%present file is the latest for the cell; simply append
                presentCell(size(presentCell,2)+1)=handles.cell;
            elseif (index==1) % present file if the first one for the cell
               newpresentCell(1,1)=handles.cell;
               newpresentCell(1,2:size(presentCell,2)+1)=presentCell;
               presentCell=newpresentCell;
            else
               newpresentCell(1,1:index-1)=presentCell(1:index-1);
               newpresentCell(1,index)=handles.cell;
               newpresentCell(1,index+1:size(presentCell,2)+1)=presentCell(index:end);
               presentCell=newpresentCell;
            end
        else %override the old cell with the new one
            presentCell(1,index)=handles.cell;        
        end   
            
    save ([handles.cellPath,handles.cellName], 'presentCell');
    end
    handles.cell_allrecs=presentCell; 
end

% --- Executes on button press in buttonLoadAnnotation.
function buttonLoadAnnotation_Callback(hObject, eventdata, handles)
% hObject    handle to buttonLoadAnnotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get the new filename.
[file,path] = uigetfile('*.mat', 'Choose an audio annotation .mat file to load:');

if isequal(file,0) || isequal(path,0)
else
    %if they don't hit cancel..
    if(isfield(handles,'annotationFileName'))
        button = questdlg('there is already an annotation loaded. You want to save this and use another?')
        if (button=='Yes')
            
            %already have an annotation object loaded, then save clear current annotation.
            saveAnnotation(handles);
            guidata(hObject, handles);
            buttonClearAll_Callback(hObject, eventdata, handles);
            handles = guidata(hObject);
            handles.audioAnnotation.delete;
            handles.audioAnnotation = mhashtable;
        end
    end
    handles.annotationFileName = [path,filesep,file];
    audioAnnotation = aaLoadHashtable(handles.annotationFileName);
    handles.audioAnnotation = audioAnnotation;
    if(handles.audioAnnotation.containsKey(handles.filename))
        currAnnot = handles.audioAnnotation.get(handles.filename);
        set(handles.popupDrug,'Value',currAnnot.drugindex);
    end
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
    
end
guidata(hObject, handles);

% --- Executes on buttonUnknown press in buttonSubsong.
function buttonSubsong_Callback(hObject, eventdata, handles)
% hObject    handle to buttonSubsong (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles = labelSyllable(handles, 100);
% guidata(hObject, handles);

% --- Executes on buttonUnknown press in buttonNoise.
function buttonNoise_Callback(hObject, eventdata, handles)
% hObject    handle to buttonNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles = labelSyllable(handles, 101);
% guidata(hObject, handles);

% --- Executes on buttonUnknown press in buttonCall.
function buttonCall_Callback(hObject, eventdata, handles)
% hObject    handle to buttonCall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles = labelSyllable(handles, 102);
% guidata(hObject, handles);

% --- Executes on buttonUnknown press in buttonUnknown.
function buttonExitCell_Callback(hObject, eventdata, handles)
% hObject    handle to buttonUnknown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles = labelSyllable(handles, 103);
% guidata(hObject, handles);
handles=saveCell(handles);
if isfield(handles,'cell_allrecs')
    handles = rmfield(handles, 'cell_allrecs');
end
guidata(hObject, handles);


% --- Executes on buttonUnknown press in buttonLoadExper.
function buttonLoadExper_Callback(hObject, eventdata, handles)
% hObject    handle to buttonLoadExper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, path] = uigetfile;
try,
    load([path,filesep,file]);
    exper.dir = path;
    handles.exper = exper;
    handles.chan = exper.sigCh(1);
    saveExper(handles.exper);
    if(length(exper.sigCh)>1)
        set(handles.popupChannel,'String',exper.sigCh);
        set(handles.popupChannel,'Value',1);
    else
        set(handles.popupChannel,'String',exper.audioCh);
        set(handles.popupChannel,'Value',1);
    end
    cla(handles.axesNavBar);
    cla(handles.axesSpecgram);
    cla(handles.axesPower);
    cla(handles.axesSignal);
catch,
end
guidata(hObject, handles);
%updateTemplates(handles);


function editFileNum_Callback(hObject, eventdata, handles)
% hObject    handle to editFileNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFileNum as text
%        str2double(get(hObject,'String')) returns contents of editFileNum as a double
filenum = get(handles.editFileNum, 'String');
filenum = str2num(filenum);
handles.startNdx = -1;
handles.endNdx = -1;
handles = loadfile(handles, filenum);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editFileNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFileNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on buttonUnknown press in buttonPrevFile.
function buttonPrevFile_Callback(hObject, eventdata, handles)
% hObject    handle to buttonPrevFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filenum = handles.filenum - 1;
handles.startNdx = -1;
handles.endNdx = -1;
handles = loadfile(handles, filenum);
guidata(hObject, handles);

% --- Executes on buttonUnknown press in buttonNextFile.
function buttonNextFile_Callback(hObject, eventdata, handles)
% hObject    handle to buttonNextFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filenum = handles.filenum + 1;
handles.startNdx = -1;
handles.endNdx = -1;
handles = loadfile(handles, filenum);
guidata(hObject, handles);


% --- Executes on buttonUnknown press in buttonDeleteSyll.
function buttonDeleteSyll_Callback(hObject, eventdata, handles)
% hObject    handle to buttonDeleteSyll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.selectedSyll ~= -1)
    nSyll = handles.selectedSyll;
    handles = deleteSyllableRectAndTxt(handles, nSyll);
    if(size(handles.txtHandles,1) >= nSyll)
        handles.txtHandles(nSyll,:) = [];
    end
    if(size(handles.rectHandles,1) >= nSyll)
        handles.rectHandles(nSyll,:) = [];
    end
    currAnnotation = handles.audioAnnotation.get(handles.filename)
    currAnnotation.segAbsStartTimes(nSyll) = [];
    currAnnotation.segFileStartTimes(nSyll) = [];
    currAnnotation.segFileEndTimes(nSyll) = [];
    currAnnotation.segType(nSyll) = [];
    selectedSyll = handles.selectedSyll;
    handles.audioAnnotation.put(handles.filename, currAnnotation);
    if(length(currAnnotation.segAbsStartTimes)== 0)
        selectedSyll = -1;
    elseif(selectedSyll > length(currAnnotation.segAbsStartTimes));
        selectedSyll = selectedSyll - 1;
    end    
    
    handles = setSelectedSyllable(handles, selectedSyll);    
    guidata(hObject, handles);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function handles = deleteSyllableRectAndTxt(handles, nSyll)
if(size(handles.txtHandles,1) >= nSyll)
    for(nTxt = 1:size(handles.txtHandles,2))
        if(ishandle(handles.txtHandles(nSyll,nTxt)))
            delete(handles.txtHandles(nSyll,nTxt));
            handles.txtHandles(nSyll,nTxt) = -1;
        end
    end
end
if(size(handles.rectHandles,1) >= nSyll)
    for(nRect = 1:size(handles.rectHandles,2))
        if(ishandle(handles.rectHandles(nSyll,nRect)))
            delete(handles.rectHandles(nSyll,nRect));
            handles.rectHandles(nSyll,nRect) = -1;
        end
    end
end
    

% --- Executes on buttonUnknown press in buttonAddSyll.
function buttonAddSyll_Callback(hObject, eventdata, handles)
% hObject    handle to buttonAddSyll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.bWaitingForAddClick = true;
handles.bWaitingForSeparationClick = false;
set(handles.buttonAddSyll, 'BackgroundColor', 'red');
guidata(hObject, handles);

% --- Executes on buttonUnknown press in buttonDeletePause.
function buttonDeletePause_Callback(hObject, eventdata, handles)
% hObject    handle to buttonDeletePause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.selectedSyll ~= -1)
    filenum = handles.filenum;
    nSyll = handles.selectedSyll;
    nSyllNext = nSyll + 1;
    currAnnotation = handles.audioAnnotation.get(handles.filename);
    numSyll = length(currAnnotation.segAbsStartTimes);
    if(nSyll<numSyll)
        %if there is a syllable beyond the selected one, then delete
        %the pause in between.
        %delete the rectangles and text for both syllables.
        handles = deleteSyllableRectAndTxt(handles,nSyll);
        handles = deleteSyllableRectAndTxt(handles,nSyllNext);
        if((size(handles.txtHandles,1)) >= nSyllNext)
            handles.txtHandles(nSyllNext,:) = [];
        end
        if((size(handles.rectHandles,1)) >= nSyllNext)
            handles.rectHandles(nSyllNext,:) = [];
        end       
        
        currAnnot = handles.audioAnnotation.get(handles.filename);
        currAnnot.segAbsStartTimes(nSyllNext) = [];
        currAnnot.segFileStartTimes(nSyllNext) = [];
        currAnnot.segFileEndTimes(nSyll) = [];
        currAnnot.segType(nSyll) = -1;
        currAnnot.segType(nSyllNext) = [];
        handles.audioAnnotation.put(handles.filename, currAnnot);
        handles = drawSyllableRect(handles, nSyll);
        if(ishandle(handles.rectHandles(nSyll,:)))
            set(handles.rectHandles(nSyll,:),'FaceColor','green','EdgeColor','green');
        end
        guidata(hObject, handles);
    end
end

% --- Executes on button press in buttonClearAll.
function buttonClearAll_Callback(hObject, eventdata, handles)
% hObject    handle to buttonClearAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=clearAll(handles);
guidata(hObject, handles);

function handles=clearAll(handles)

for(nSyll = 1:max(size(handles.txtHandles,1),size(handles.rectHandles,1)))
    handles = deleteSyllableRectAndTxt(handles, nSyll);
end
handles.txtHandles = [];
handles.rectHandles = [];
handles.audioAnnotation.remove(handles.filename);
handles.selectedSyll = -1;
%guidata(hObject, handles);


% --- Executes on button press in buttonAutoSeg.
function buttonAutoSeg_Callback(hObject, eventdata, handles)
% hObject    handle to buttonAutoSeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%if there has been segmentation before delete all boxes and numbers
if(size(handles.rectHandles,1)>1)
    handles=clearAll(handles);
end
fs = handles.exper.desiredInSampRate;
startTime = (handles.startNdx-1)/handles.fs;
startNdx = handles.startNdx;
endNdx = handles.endNdx;
audio_select=handles.audio(startNdx:endNdx);
[syllStartTimes, syllEndTimes, noiseEst, noiseStd, soundEst, thresEdge,thresSyll, soundStd] = aSAP_segSyllablesFromRawAudioBence(audio_select, fs, handles.edges, handles.threshold);
syllStartTimes=syllStartTimes+startTime;
syllEndTimes=syllEndTimes+startTime;
handles.thresSyll=thresSyll;
handles.thresEdge=thresEdge;
set(handles.editSyllAbs, 'String',num2str(thresSyll,4));
set(handles.editEdgeAbs, 'String',num2str(thresEdge,4));

time = extractExperFilenumTime(handles.exper,handles.filenum);
segType=repmat(-1,length(syllStartTimes), 1);
jitterMS = .010;
matchThreshold = 0.025;
%if there is a template, compare each identified syllable to the template
%syllables

if (isfield(handles,'templateFeatures') && get(handles.checkboxAutoID, 'Value'))
    matchScores(1:size(handles.templateFeatures,2))=0;
    for nSyll=1:length(syllStartTimes)
        audio=handles.audio(round(((syllStartTimes(nSyll)-0.015)*fs)):round(((syllEndTimes(nSyll)+0.005)*fs)));
        [SAPFeats, m_spec_deriv]  = aSAP_generateASAPFeatures(audio, fs);
        for nTemp = 1:size(handles.templateFeatures,2)
            if (~isempty(handles.templateFeatures(nTemp).lengthAudio))
                matchScores(nTemp) = aSAP_computeMatchScore1(SAPFeats, handles.templateFeatures(nTemp), jitterMS, false);           
            end
        end
        matchScores
        [topScore,topScoreIndex]=max(matchScores);
        if (topScore>matchThreshold)
            segType(nSyll)=topScoreIndex;
        end
    end
end
currAnnot.exper = handles.exper;
currAnnot.filenum = handles.filenum;
currAnnot.segAbsStartTimes = time + (syllStartTimes/(24*60*60));
currAnnot.segFileStartTimes = syllStartTimes;
currAnnot.segFileEndTimes = syllEndTimes;
currAnnot.segType = segType;
currAnnot.fs = fs;
currAnnot.drugstatus=handles.drugstatus;
currAnnot.syllthresh=handles.threshold;
currAnnot.edgethresh=handles.edges;
currAnnot.abssyllthresh=handles.thresSyll;
currAnnot.absedgethresh=handles.thresEdge;
currAnnot.drugindex=handles.drugindex;
currAnnot.drugstatus=handles.drugstatus;
handles.audioAnnotation.put(handles.filename, currAnnot);

%draw syllable information
handles.rectHandles = [];
handles.txtHandles = [];
handles = drawNavBarSyllableRects(handles);
handles = drawZoomSyllableRects(handles);
guidata(hObject, handles);


% --- Executes on selection change in popupChannel.
function popupChannel_Callback(hObject, eventdata, handles)
% hObject    handle to popupChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupChannel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupChannel
channels = get(handles.popupChannel, 'String');
nValue = get(handles.popupChannel, 'Value');
handles.chan =channels(nValue);
handles.cell_no=1;
axes(handles.axesSignal);
try
    handles.sig = loadData(handles.exper, handles.filenum, handles.chan);
    plotTimeSeriesQuick(handles.time(handles.startNdx:handles.endNdx), handles.sig(handles.startNdx:handles.endNdx));
    set(handles.axesSignal, 'ButtonDownFcn', '');
    axes tight;
catch
    gca;
end
handles=ShowSpikes(handles);
guidata(hObject, handles);

function chan = getChan(handles)
channels = get(handles.popupChannel, 'String');
nValue = get(handles.popupChannel, 'Value');
chan = channels(nValue);

% --- Executes during object creation, after setting all properties.
function popupChannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in buttonSaveTemplates.
function buttonSaveTemplates_Callback(hObject, eventdata, handles)
% hObject    handle to buttonSaveTemplates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(isfield(handles,'templates'))
    templates=handles.templates;
    templateFeatures=handles.templateFeatures;
    [FileName,PathName]=uiputfile('*.mat','Save Template As:','template');
    save([PathName FileName],'templates', 'templateFeatures');
else
    warndlg('No templates to save.');
    uiwait;
end

% --- Executes on button press in buttonLoadTemplates.
function buttonLoadTemplates_Callback(hObject, eventdata, handles)
% hObject    handle to buttonLoadTemplates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiload;
handles.templates = templates;
handles.templateFeatures = templateFeatures;
updateTemplates(handles);
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = labelSyllable(handles, label)
if(handles.filenum > 0 && handles.selectedSyll > 0)
    filenum = handles.filenum;
    nSyll = handles.selectedSyll;
    fs = handles.fs;
    
    %set the label value
    currAnnot = handles.audioAnnotation.get(handles.filename);
    currAnnot.segType(nSyll) = label;
    handles.audioAnnotation.put(handles.filename, currAnnot);
    
    %delete any old text handles.
    if((size(handles.txtHandles,1)) >= nSyll)
        for(nTxt = 1:size(handles.txtHandles,2))
            if ishandle(handles.txtHandles(nSyll,nTxt))
                delete(handles.txtHandles(nSyll,nTxt));
                handles.txtHandles(nSyll,nTxt) = -1;
            end
        end
    end
    
    %Compute middle of syllable
    midSyll = (currAnnot.segFileEndTimes(nSyll) + currAnnot.segFileStartTimes(nSyll)) / 2;
    
    %create a new text
    axes(handles.axesSpecgram);
    handles.txtHandles(nSyll, 1) = text(midSyll,mean(ylim),num2str(label));
    axes(handles.axesPower);
    handles.txtHandles(nSyll, 2) = text(midSyll,mean(ylim),num2str(label));
    axes(handles.axesSignal);
    handles.txtHandles(nSyll, 3) = text(midSyll,mean(ylim),num2str(label));
    set(handles.txtHandles(nSyll,1:3), 'Color', 'black');
    set(handles.txtHandles(nSyll,1:3), 'FontSize', 14);
    set(handles.txtHandles(nSyll,1:3), 'FontWeight', 'bold');
    set(handles.txtHandles(nSyll,1:3), 'HitTest', 'off');
    set(handles.txtHandles(nSyll,1:3), 'HorizontalAlignment', 'center'); 
 
    %change rect colors to blue
    if(ishandle(handles.rectHandles(nSyll,:)))
        set(handles.rectHandles(nSyll,:),'FaceColor','blue','EdgeColor','none');
        set(handles.rectHandles(nSyll,2),'EdgeColor','blue');
    end
    
    %set the current lable text
    set(handles.textCurrentLabel,'String', ['Current Label: ', num2str(label)]);
    %select next syllable and shift specgram axis.
    handles = setSelectedSyllable(handles, handles.selectedSyll + 1);
    if(currAnnot.segFileEndTimes(nSyll) > (handles.endNdx-1)/fs)
        width = handles.endNdx - handles.startNdx;
        handles.startNdx = max(0,round(round((midSyll*fs)+1) - width/2));
        handles.endNdx = min(length(handles.audio), round(round((midSyll*fs)+1) + width/2));
        handles = updateZoomOrPan(handles);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = drawSyllableRect(handles, nSyll, bInsert)
fs = handles.fs;
windowStartTime = (handles.startNdx-1)/fs;
windowEndTime = (handles.endNdx-1)/fs;
currAnnot = handles.audioAnnotation.get(handles.filename);
startTime = currAnnot.segFileStartTimes(nSyll);
endTime = currAnnot.segFileEndTimes(nSyll);
label = currAnnot.segType(nSyll);
if(label ~= -1)
    color = [0 0 1];
else
    color = [1 0 0];
end
if(handles.selectedSyll == nSyll)
    color = [0 1 0];
end

if(~exist('bInsert') || ~bInsert)
    handles = deleteSyllableRectAndTxt(handles, nSyll);
end

x = [startTime, startTime, endTime, endTime, startTime];
axes(handles.axesNavBar)
lims = ylim;
y = [lims(1), lims(2), lims(2), lims(1), lims(1)];

if(exist('bInsert') && bInsert)
    handles.rectHandles = [handles.rectHandles([1:nSyll-1],:);repmat(-1,1,size(handles.rectHandles,2));handles.rectHandles([nSyll:end],:)];
    handles.txtHandles = [handles.txtHandles([1:nSyll-1],:);repmat(-1,1,size(handles.txtHandles,2));handles.txtHandles([nSyll:end],:)];
end

handles.rectHandles(nSyll, 1) = patch(x,y, color, 'EdgeColor', 'none', 'FaceAlpha', .2);
set(handles.rectHandles(nSyll,1), 'HitTest', 'off'); 

if(startTime<=windowEndTime & endTime>=windowStartTime)    
    axes(handles.axesSpecgram)
    lims = ylim;
    y = [lims(1), lims(2), lims(2), lims(1), lims(1)];
    handles.rectHandles(nSyll, 2) = patch(x,y, color, 'EdgeColor', color, 'FaceAlpha', 0);

    axes(handles.axesPower)
    lims = ylim;
    y = [lims(1), lims(2), lims(2), lims(1), lims(1)];
    handles.rectHandles(nSyll, 3) = patch(x,y, color, 'EdgeColor', 'none', 'FaceAlpha', .2);

    axes(handles.axesSignal)
    lims = ylim;
    y = [lims(1), lims(2), lims(2), lims(1), lims(1)];
    handles.rectHandles(nSyll, 4) = patch(x,y, color, 'EdgeColor', 'none', 'FaceAlpha', .2);
    
    set(handles.rectHandles(nSyll,2:4), 'HitTest', 'off'); 

    if(label ~= -1) 
        midSyll = mean([startTime,endTime]);
        %create a new text handles
        axes(handles.axesSpecgram);
        handles.txtHandles(nSyll, 1) = text(midSyll,mean(ylim),num2str(label));
        axes(handles.axesPower);
        handles.txtHandles(nSyll, 2) = text(midSyll,mean(ylim),num2str(label));
        axes(handles.axesSignal);
        handles.txtHandles(nSyll, 3) = text(midSyll,mean(ylim),num2str(label));
        set(handles.txtHandles(nSyll,:), 'Color', 'black');
        set(handles.txtHandles(nSyll,:), 'FontSize', 14);
        set(handles.txtHandles(nSyll,:), 'FontWeight', 'bold');
        set(handles.txtHandles(nSyll,:), 'HitTest', 'off');   
        set(handles.txtHandles(nSyll,:), 'HorizontalAlignment', 'center');         
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = loadfile(handles, filenum)
saveAnnotation(handles);
handles=saveCell(handles);
handles.sorted_spikes=[];
handles.unsorted_spikes=[];
%handles.cell.spont=[];
if(filenum == handles.filenum)
    return;
end
try
    fs = handles.exper.desiredInSampRate;
    handles.filename = getExperDatafile(handles.exper,filenum,handles.exper.audioCh);
    handles.audio = loadAudio(handles.exper, filenum);
    handles.sig = loadData(handles.exper, filenum, handles.chan);
catch
    warndlg('File number specified does not exist.');
    set(handles.editFileNum,'String',num2str(handles.filenum));
    uiwait;
    return;
end
handles.filenum = filenum;
set(handles.editFileNum,'String',num2str(handles.filenum));
handles.time = linspace(0, (length(handles.audio)-1)/fs, length(handles.audio));
if(handles.endNdx - handles.startNdx > .01)
    width = min(handles.endNdx - handles.startNdx, length(handles.audio));
else
    width = length(handles.audio);
end
handles.startNdx = 1;
handles.endNdx = width;
handles.fs = fs;
[pow, filtAud] = aSAP_getLogPower(handles.audio, fs);
handles.audio = filtAud;
handles.power = pow;
handles.selectedSyll = -1;
handles.rectHandles = [];
handles.txtHandles = [];
axes(handles.axesNavBar);
plotTimeSeriesQuick(handles.time, handles.audio);
%set(handles.axesNavBar,'Color','black');

%Draw rectangles on the nav bar...
if(handles.audioAnnotation.containsKey(handles.filename))
    currAnnot = handles.audioAnnotation.get(handles.filename);
    handles.drugstatus=currAnnot.drugstatus;
    handles.threshold=currAnnot.syllthresh;
    handles.edges=currAnnot.edgethresh;
    handles.thresSyll=currAnnot.abssyllthresh;
    handles.thresEdge=currAnnot.absedgethresh;
    handles.drugindex=currAnnot.drugindex;
    set(handles.editEdgeThreshold,'String',num2str(handles.edges));
    set(handles.editPowerThreshold,'String',num2str(handles.threshold));
    set(handles.editSyllAbs, 'String',num2str(handles.thresSyll));
    set(handles.editEdgeAbs, 'String',num2str(handles.thresEdge));
    set(handles.popupDrug,'Value',handles.drugindex);
    contents = get(handles.popupDrug,'String');% returns popupDrug contents as cell array
    if (~strcmp(currAnnot.drugstatus,contents{handles.drugindex}))
        errordlg('drug definitions have changed')
    end
 
    syllStartTimes = currAnnot.segFileStartTimes;
    syllEndTimes = currAnnot.segFileEndTimes;
    x = [syllStartTimes; syllStartTimes; syllEndTimes; syllEndTimes; syllStartTimes];
    axes(handles.axesNavBar)
    lims = ylim;
    y = [lims(1); lims(2); lims(2); lims(1); lims(1)];
    for nSyll = 1:length(syllStartTimes)
        if(currAnnot.segType(nSyll) == -1)
            color = 'red';
        else
            color = 'blue';
        end
        handles.rectHandles(nSyll, 1) = patch(x(:,nSyll),y, color, 'EdgeColor', 'none', 'FaceAlpha', .2);
    end
    set(handles.rectHandles(:,1), 'HitTest', 'off');
end
%fill up cell

handles.cell.birdname=handles.exper.birdname;
handles.cell.channel=handles.chan;
handles.cell.cell_no=handles.cell_no;
handles.cell.filenum=handles.filenum;
handles.cell.filename=handles.filename;
handles = updateZoomOrPan(handles);
handles=ShowSpikes(handles);

% 
% axes(handles.axesSpecgram);
% displaySpecgramQuick(handles.audio, fs, [0,10000]);
% 
% axes(handles.axesPower);
% plotTimeSeriesQuick(handles.time, handles.power);
% set(handles.axesPower, 'ButtonDownFcn', '');
% 
% axes(handles.axesSignal);
% plotTimeSeriesQuick(handles.time, handles.sig);
% set(handles.axesSignal,'Color','black');
% set(handles.axesSignal, 'ButtonDownFcn', '');
% 
% linkaxes([handles.axesSignal, handles.axesPower, handles.axesSpecgram],'x');
% set(get(handles.axesNavBar,'Parent'), 'KeyPressFcn', @cb_keypress);
% set(handles.axesNavBar, 'ButtonDownFcn', @cb_navbar_click);
% set(handles.axesSpecgram, 'ButtonDownFcn', @cb_specgram_click);
% 
% 
% if((filenum <= length(handles.audioAnnotation)) && ~isequal(handles.audioAnnotation{filenum},[]))
%     for(nSyll = 1:length(handles.audioAnnotation{filenum}.segFileStartTimes))
%         handles = drawSyllableRect(handles, nSyll);
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function handles = setSelectedSyllable(handles, newSelectedSyll)
%make sure selection is in bounds
currAnnot = handles.audioAnnotation.get(handles.filename);
if(newSelectedSyll == -1)
    handles.selectedSyll = -1;
    return;
elseif(newSelectedSyll < 1 | newSelectedSyll > length(currAnnot.segAbsStartTimes))
    return
end

if(handles.selectedSyll ~= -1 && handles.selectedSyll <= length(currAnnot.segAbsStartTimes))
    %unhighlight the old selected syll
    if(currAnnot.segType(handles.selectedSyll) == -1)
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
filenum = handles.filenum;
fs = handles.exper.desiredInSampRate;
startTime = (handles.startNdx-1)/handles.fs;
endTime = (handles.endNdx-1)/handles.fs;
syllEndTime = currAnnot.segFileEndTimes(selectedSyll);
syllStartTime = currAnnot.segFileStartTimes(selectedSyll);
if(syllEndTime > endTime || syllStartTime < startTime)
    width = handles.endNdx - handles.startNdx;
    mid = (syllStartTime*fs) + 1;
    handles.startNdx = round(mid - width/2);
    handles.endNdx = round(mid + width/2);
    if(handles.startNdx<1)
        handles.startNdx = 1;
        handles.endNdx = width + handles.startNdx;
    end
    if(handles.endNdx > length(handles.audio))
        handles.startNdx = length(handles.audio) - width;
        handles.endNdx = length(handles.audio);
    end        
    handles = updateZoomOrPan(handles);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateTemplates(handles)
cla(handles.axesTemplates);

if(isfield(handles, 'templates'))
    stdfs = 44100;
    allwavs = [];
    edges = [0];
    for(nWav = 1:length(handles.templates.wavs))
        fs = handles.templates.wavs(nWav).fs;
        wav = handles.templates.wavs(nWav).wav;
        label = handles.templates.wavs(nWav).segType;
        wav = resample(wav, stdfs, fs);
        allwavs = [allwavs; wav];
        edges = [edges, edges(end) + (length(wav) - 1)/stdfs];
    end
    
    %zero pad to 1 sec in length.
    if(edges(end) < 1)
        allwavs(end+1:stdfs) = 0;
    end
    
    %display the templates
    axes(handles.axesTemplates);
    displaySpecgramQuick(allwavs, stdfs, [0,10000]);
    
    %[SAPFeats, m_spec_deriv] = aSAP_generateASAPFeatures(allwavs, fs, Parameters);
    %aSAP_displaySpectralDerivative(m_spec_deriv, ...
    %    Parameters, 0, 0, -inf, false, ...
    %    0, 1, .11, ...
    %    false, false);
    
    %seperate with lines and label.
    for(nWav = 1:length(edges)-1)
        l = line([edges(nWav),edges(nWav)], ylim);
        set(l,'Color','k');
        set(l,'LineWidth', 2);
        mid = mean(edges([nWav,nWav+1]));
        if(handles.templates.wavs(nWav).segType ~= -1)
            t = text(mid, mean(ylim), num2str(handles.templates.wavs(nWav).segType));
        end
    end
    
    set(handles.axesTemplates, 'ButtonDownFcn', @cb_templates_click);
end

function handles=updateTemplateFeatures(handles,nselect,action)

%nselect is the syllable number (from handles.templates) that is being added or deleted, 
%action refers to whether a syllable should be added or deleted from the template
fs=40000;
if(isfield(handles, 'templates'))
    if(nselect==-1)
        return;
    end
    %add the featues to the feature template
    audio = handles.templates.wavs(nselect).wav;
    label = handles.templates.wavs(nselect).segType;
    if(label>0 && label<11)
        if (strcmp(action,'add'))
        %if a syllable calculate features
               
            [SAPFeats, m_spec_deriv]  = aSAP_generateASAPFeatures(audio, fs);
            handles.templateFeatures(label)=SAPFeats;
%             segTypes=getFieldVector(handles.templates.wavs,segType);
%             numberSeg=find(segTypes==label);
%             if (numberSeg==1)
%                 handles.templateFeatures(label)=SAPFeats;
%             else
%                 handles=addSyllToTemplate(handles,label,numberSeg, SAPFeats);
%             end
        
        elseif (strcmp(action,'delete'))
            label = handles.templates.wavs(nselect).segType;
            handles.templateFeatures(label)=[];
        end
    end
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = updateZoomOrPan(handles)

fs = handles.exper.desiredInSampRate;
startTime = (handles.startNdx-1)/handles.fs;
endTime = (handles.endNdx-1)/handles.fs;
startNdx = handles.startNdx;
endNdx = handles.endNdx;

%fix navigator rect
axes(handles.axesNavBar);
if(ishandle(handles.navRect))
    delete(handles.navRect);
end
x = [startTime, startTime, endTime, endTime, startTime];
y = ylim;
y = [y(1), y(2), y(2), y(1), y(1)];
handles.navRect = line(x,y);
set(handles.navRect,'Color', 'red');
set(handles.navRect,'LineWidth', 2);
set(handles.navRect,'HitTest','off');

SpecAudio = handles.audio*1;

%replot everything
axes(handles.axesSpecgram);
%displaySpecgramQuick(handles.audio(startNdx:endNdx), fs,[0,10000],[],startTime);
displaySpecgramQuick(SpecAudio(startNdx:endNdx), fs,[0,10000],[],startTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Below is the patched in code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% colormap(bbyrp(handles))
% handles.spectrumplots = [];
% handles.bgcolor = [.5 .5 .5];
% [p freq t] = power_spectrum(handles.audio(startNdx:endNdx),fs,1);
% handles.satur = 25; %looks like this number was always 25
% set(gca,'color',handles.bgcolor);
% %subplot(handles.axesSpecgram);
% %h = imagesc(t+lst(c)/handles.fs,freq/1000,atan(p/(handles.satur*10^10)));
% h = imagesc(t,freq/1000,atan(p/(handles.satur*10^10)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xlabel('');
ylabel('');
set(handles.axesSpecgram,'XTick',[]);
axes(handles.axesPower);
plotTimeSeriesQuick(handles.time(startNdx:endNdx), handles.power(startNdx:endNdx));
set(handles.axesPower,'ButtonDownFcn','');
set(handles.axesPower,'XTick',[]);

axes(handles.axesSignal);
plotTimeSeriesQuick(handles.time(startNdx:endNdx), handles.sig(startNdx:endNdx));
%set(handles.axesSignal,'Color','black');
set(handles.axesSignal,'ButtonDownFcn','');

linkaxes([handles.axesSignal, handles.axesPower, handles.axesSpecgram, handles.axesRasters],'x');
set(get(handles.axesNavBar,'Parent'), 'KeyPressFcn', @cb_keypress);
set(handles.axesNavBar, 'ButtonDownFcn', @cb_navbar_click);
set(handles.axesSpecgram, 'ButtonDownFcn', @cb_specgram_click);
set(handles.axesRasters, 'ButtonDownFcn', @cb_rasters_click);

handles.rectHandles(:,2:4) = -1;
handles.txtHandles(:,1:3) = -1;
handles = drawZoomSyllableRects(handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
 function cb_keypress(hObject, evnt)
% hObject    handle to buttonAutoSeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

%get modifier booleans
bShift = false;
bControl = false;
if(length(evnt.Modifier) >= 1)
    for(nMod = 1:length(evnt.Modifier))
        bShift = bShift || strcmp(evnt.Modifier{nMod}, 'shift');
        bControl = bControl || strcmp(evnt.Modifier{nMod}, 'control');
    end
end

[label, ok] = str2num(evnt.Key);
if(ok)
    if(bShift)
        label = label + 200;
    elseif(bControl)
        label = label + 300
    end
    handles = labelSyllable(handles, label);
else
    if(strcmp(evnt.Key,'u')) 
        handles = labelSyllable(handles, 103);
    elseif(strcmp(evnt.Key,'s')) 
        handles = labelSyllable(handles, 100);
    elseif(strcmp(evnt.Key,'n')) 
        handles = labelSyllable(handles, 101);
    elseif(strcmp(evnt.Key,'c')) 
        handles = labelSyllable(handles, 102);
    elseif(strcmp(evnt.Key,'period')) 
        width = handles.endNdx - handles.startNdx;
        if(handles.endNdx + width/2 > length(handles.audio))
            handles.startNdx = length(handles.audio) - width;
            handles.endNdx = length(handles.audio);
        else
            handles.startNdx = round(handles.startNdx + width/1.5);
            handles.endNdx = round(handles.endNdx + width/1.5);
        end
        handles = updateZoomOrPan(handles);
    elseif(strcmp(evnt.Key,'comma'))
        width = handles.endNdx - handles.startNdx;
        if(handles.startNdx - width/2 < 1)
            handles.startNdx = 1;
            handles.endNdx = width;
        else
            handles.startNdx = round(handles.startNdx - width/2);
            handles.endNdx = round(handles.endNdx - width/2);
        end
        handles = updateZoomOrPan(handles);
    elseif(strcmp(evnt.Key,'space')) 
        if(bShift)
            handles = setSelectedSyllable(handles, handles.selectedSyll - 1)   
        else
            handles = setSelectedSyllable(handles, handles.selectedSyll + 1)
        end
    elseif(strcmp(evnt.Key,'p')) 
        soundsc(handles.audio(handles.startNdx:handles.endNdx), handles.exper.desiredInSampRate);
    elseif(strcmp(evnt.Key,'delete'))
        guidata(hObject, handles);
        buttonDeleteSyll_Callback(hObject, [], handles)
        handles = guidata(hObject);
    elseif(strcmp(evnt.Key,'m'))
        guidata(hObject, handles);
        buttonDeletePause_Callback(hObject, [], handles)
        handles = guidata(hObject);
    end
end
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function cb_navbar_click(hObject, evnt)
handles = guidata(hObject);
handles.fs = handles.exper.desiredInSampRate;
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
        if(handles.endNdx > length(handles.audio))
            handles.startNdx = length(handles.audio) - width;
            handles.endNdx = length(handles.audio);
        end        
        handles = updateZoomOrPan(handles);
    else
        handles.startNdx = max(1,round((p1(1)*fs) + 1));
        handles.endNdx = min(length(handles.audio), round(((p1(1) + offset(1))*fs) + 1));
        handles = updateZoomOrPan(handles);
    end
end
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_specgram_click(hObject, evnt)
handles = guidata(hObject);
fs = handles.exper.desiredInSampRate;
filenum = handles.filenum;
mouseMode = get(get(hObject,'Parent'), 'SelectionType');
clickLocation = get(handles.axesSpecgram, 'CurrentPoint');
axes(handles.axesSpecgram);

if(handles.bWaitingForAddClick)
    %handle add syllable click...
    handles.bWaitingForAddClick = false;
    set(handles.buttonAddSyll, 'BackgroundColor', 'white');
    
    %create audio annotation if it doesn't already exist:
    if(~handles.audioAnnotation.containsKey(handles.filename))
        currAnnot.exper = handles.exper;
        currAnnot.filenum = handles.filenum;
        currAnnot.segAbsStartTimes = [];
        currAnnot.segFileStartTimes = [];
        currAnnot.segFileEndTimes = [];
        currAnnot.segType = [];
        currAnnot.fs = handles.fs;
        currAnnot.drugstatus=handles.drugstatus;
        currAnnot.syllthresh=handles.threshold;
        currAnnot.edgethresh=handles.edges;
        currAnnot.abssyllthresh=handles.thresSyll;
        currAnnot.absedgethresh=handles.thresEdge;
        currAnnot.drugindex=handles.drugindex;
        currAnnot.drugstatus=handles.drugstatus;
        
    else
        currAnnot = handles.audioAnnotation.get(handles.filename);
    end

    %get region for addition
    rect = rbbox;
    endPoint = get(gca,'CurrentPoint'); 
    point1 = clickLocation(1,1:2);              % extract x and y
    point2 = endPoint(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    offset = abs(point1-point2);         % and dimensions
    
    if(get(handles.checkboxEdgefind, 'Value'))
        %adjust p1 and offset accordingly.
        startNdx = round((p1(1)*handles.fs) + 1);
        endNdx = round(((p1(1)+ offset(1))*handles.fs) + 1);
        edgeThres = get(handles.editEdgeAbs, 'String');
        [edgeThres,ok] = str2num(edgeThres);
        if(ok)
            ndx = find(handles.power(startNdx:endNdx) > edgeThres);
            if(length(ndx)>=2)
                p1(1) = ((startNdx + ndx(1) - 1) - 1) / handles.fs;
                offset(1) = (((startNdx + ndx(end) - 1) - 1) / handles.fs) - p1(1);
            end
        end
    end
    
    %be sure the region doesn't overlap another syllable.
    syllStartTimes = currAnnot.segFileStartTimes;
    syllEndTimes = currAnnot.segFileEndTimes;
    nSyllAfter = find(p1(1) + offset(1) < [syllStartTimes,Inf]); nSyllAfter = nSyllAfter(1);
    nSyllBefore = find(p1(1) > [-Inf, syllEndTimes]); nSyllBefore = nSyllBefore(end)-1;
    if(nSyllBefore - nSyllAfter == -1)
        time = extractExperFilenumTime(handles.exper,handles.filenum);
        currAnnot.segAbsStartTimes = [currAnnot.segAbsStartTimes([1:nSyllBefore]),time + (p1(1)/(24*60*60)),currAnnot.segAbsStartTimes([nSyllAfter:end])];
        currAnnot.segFileStartTimes = [currAnnot.segFileStartTimes([1:nSyllBefore]),p1(1),currAnnot.segFileStartTimes([nSyllAfter:end])];
        currAnnot.segFileEndTimes = [currAnnot.segFileEndTimes([1:nSyllBefore]),p1(1) + offset(1),currAnnot.segFileEndTimes([nSyllAfter:end])];
        currAnnot.segType = [currAnnot.segType([1:nSyllBefore]);-1;currAnnot.segType([nSyllAfter:end])];
        handles.audioAnnotation.put(handles.filename, currAnnot);
        handles = drawSyllableRect(handles, nSyllBefore + 1, true);

    else
        warndlg('Specified segment overlaps another segment.');
        uiwait;
    end    
elseif(handles.bWaitingForSeparationClick)
    handles.bWaitingForSeparationClick = false;
elseif(handles.audioAnnotation.containsKey(handles.filename))
    %Otherwise check for syllable to select.
    currAnnot = handles.audioAnnotation.get(handles.filename);
    nSelect = find(clickLocation(1,1) >= currAnnot.segFileStartTimes & ...
                   clickLocation(1,1) <= currAnnot.segFileEndTimes);
    if(length(nSelect) == 1)
        if(nSelect ~= handles.selectedSyll)
            handles = setSelectedSyllable(handles, nSelect);
        end
        if(strcmp(mouseMode, 'open'))
            if(~isfield(handles, 'templates'))
                handles.templates = [];
                handles.templates.wavs = [];
                %handles.templateFeatures=[];
            else
                    if (size(handles.templates.wavs,2)~=0)
                        segTypes=getFieldVector(handles.templates.wavs,'segType');
                        if (~isempty(find(segTypes==currAnnot.segType(handles.selectedSyll))))
                            warndlg('Delete the existing syllable with the same Annotation first');
                        return;
                        end
                    end
            end


            startTime = currAnnot.segFileStartTimes(handles.selectedSyll);
            endTime = currAnnot.segFileEndTimes(handles.selectedSyll);
            if (isempty(handles.templates.wavs))
                handles.templates.wavs(1).filename = getExperDatafile(handles.exper, handles.filenum, handles.chan);
                handles.templates.wavs(1).startTime = startTime;
                handles.templates.wavs(1).endTime = endTime;
                handles.templates.wavs(1).fs = fs;
                handles.templates.wavs(1).wav = handles.audio(round((startTime*fs)+1):round((endTime*fs)+1));
                handles.templates.wavs(1).segType = currAnnot.segType(handles.selectedSyll); 
            else
                segType=getFieldVector(handles.templates.wavs,'segType');
                index=find(segType==currAnnot.segType(handles.selectedSyll));
                if (isempty(index))
                    index=find(segType<currAnnot.segType(handles.selectedSyll),1,'last')+1;
                end
                if (isempty(index))
                    index=1;
                end
                for i=length(handles.templates.wavs)+1:-1:index+1
                    handles.templates.wavs(i).filename = handles.templates.wavs(i-1).filename;
                    handles.templates.wavs(i).startTime = handles.templates.wavs(i-1).startTime;
                    handles.templates.wavs(i).endTime = handles.templates.wavs(i-1).endTime;
                    handles.templates.wavs(i).fs = handles.templates.wavs(i-1).fs;
                    handles.templates.wavs(i).wav = handles.templates.wavs(i-1).wav 
                    handles.templates.wavs(i).segType = handles.templates.wavs(i-1).segType;
                end
                handles.templates.wavs(index).filename = getExperDatafile(handles.exper, handles.filenum, handles.chan);
                handles.templates.wavs(index).startTime = startTime;
                handles.templates.wavs(index).endTime = endTime;
                handles.templates.wavs(index).fs = fs;
                handles.templates.wavs(index).wav = handles.audio(round((startTime*fs)+1):round((endTime*fs)+1));
                handles.templates.wavs(index).segType = currAnnot.segType(handles.selectedSyll); 
            end
            
            updateTemplates(handles);
            handles=updateTemplateFeatures(handles,size(handles.templates.wavs,2),'add');
        
            
        end
    end
end
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_templates_click(hObject, evnt)
handles = guidata(hObject);
fs = handles.exper.desiredInSampRate;
filenum = handles.filenum;
mouseMode = get(get(hObject,'Parent'), 'SelectionType');
clickLocation = get(handles.axesTemplates, 'CurrentPoint');
axes(handles.axesTemplates);
if(strcmp(mouseMode, 'open'))
    if(isfield(handles, 'templates'))
        stdfs = 44100;
        edges = [0];
        allwavs = [];
        for(nWav = 1:length(handles.templates.wavs))
            fs = handles.templates.wavs(nWav).fs;
            wav = handles.templates.wavs(nWav).wav;
            label = handles.templates.wavs(nWav).segType;
            wav = resample(wav, stdfs, fs);
            allwavs = [allwavs; wav];
            edges = [edges, edges(end) + (length(wav) - 1)/stdfs];
        end
        nSelect = find(clickLocation(1,1) < edges);
        if(length(nSelect) > 0)
            nSelect = nSelect(1) - 1;
            handles=updateTemplateFeatures(handles,nSelect,'delete');
            handles.templates.wavs(nSelect) = [];
            updateTemplates(handles);
            
        end
    end
end
guidata(hObject, handles);

function cb_rasters_click(hObject, evnt)
handles = guidata(hObject);
click = get(gca,'CurrentPoint');
% if (~isempty(handles.spikemarker))
%     delete(handles.spikemarker);
% end
handles.spiketime=click(1,1);
hold on
handles.spikemarker = plot(handles.spiketime, -0.1,'g:*'); 
hold off
guidata(hObject, handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = drawNavBarSyllableRects(handles)
filenum = handles.filenum;
fs = handles.exper.desiredInSampRate;
startTime = (handles.startNdx-1)/handles.fs;
endTime = (handles.endNdx-1)/handles.fs;
startNdx = handles.startNdx;
endNdx = handles.endNdx;

if(handles.audioAnnotation.containsKey(handles.filename))
    currAnnot = handles.audioAnnotation.get(handles.filename);
    syllStartTimes = currAnnot.segFileStartTimes;
    syllEndTimes = currAnnot.segFileEndTimes;
    syllType = currAnnot.segType;
    x = [syllStartTimes; syllStartTimes; syllEndTimes; syllEndTimes; syllStartTimes];   
    color = zeros(3,length(syllStartTimes));
    color(3,syllType~=-1) = 1;
    color(1,syllType==-1) = 1;
    if(handles.selectedSyll~=-1)
        color(:,handles.selectedSyll) = [0,1,0];
    end
    
    %draw nav bar rectangles.
    axes(handles.axesNavBar)
    lims = ylim;
    y = [lims(1); lims(2); lims(2); lims(1); lims(1)];
    for(nSyll = 1:length(syllStartTimes))
        handles.rectHandles(nSyll, 1) = patch(x(:,nSyll),y, color(:,nSyll)', 'EdgeColor', 'none', 'FaceAlpha', .2);
    end
    set(handles.rectHandles(:,1), 'HitTest', 'off');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = drawZoomSyllableRects(handles)
filenum = handles.filenum;
fs = handles.exper.desiredInSampRate;
startTime = (handles.startNdx-1)/handles.fs;
endTime = (handles.endNdx-1)/handles.fs;
startNdx = handles.startNdx;
endNdx = handles.endNdx;

if(handles.audioAnnotation.containsKey(handles.filename))
    currAnnot = handles.audioAnnotation.get(handles.filename);
    syllStartTimes = currAnnot.segFileStartTimes;
    syllEndTimes = currAnnot.segFileEndTimes;
    midSyll = (syllStartTimes + syllEndTimes) / 2;
    syllType = currAnnot.segType;
    x = [syllStartTimes; syllStartTimes; syllEndTimes; syllEndTimes; syllStartTimes];   
    color = zeros(3,length(syllStartTimes));
    color(3,syllType~=-1) = 1;
    color(1,syllType==-1) = 1;
    if(handles.selectedSyll~=-1)
        color(:,handles.selectedSyll) = [0,1,0];
    end
 
    %Only the visible rectangles need be drawn in the other axes.
    ndx = find(syllStartTimes < endTime & syllEndTimes> startTime);

    axes(handles.axesSpecgram)
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
    
    axes(handles.axesSignal)
    lims = ylim;
    y = [lims(1); lims(2); lims(2); lims(1); lims(1)];
    for(nSyll = ndx)
        handles.rectHandles(nSyll, 4) = patch(x(:,nSyll),y, color(:,nSyll)', 'EdgeColor', 'none', 'FaceAlpha', .2);
        if(syllType(nSyll) ~= -1)
            handles.txtHandles(nSyll, 3) = text(midSyll(nSyll),mean(ylim),num2str(syllType(nSyll)));
            set(handles.txtHandles(nSyll,3), 'Color', 'black');
            set(handles.txtHandles(nSyll,3), 'FontSize', 14);
            set(handles.txtHandles(nSyll,3), 'FontWeight', 'bold');
            set(handles.txtHandles(nSyll,3), 'HitTest', 'off'); 
            set(handles.txtHandles(nSyll,3), 'HorizontalAlignment', 'center'); 
        else
            handles.txtHandles(nSyll,3) = -1;
        end
    end
    
    set(handles.rectHandles(ndx,2:4), 'HitTest', 'off'); 
    handles.rectHandles(handles.rectHandles==0) = -1;
    handles.rectHandles(handles.txtHandles==0) = -1;
end

% --- Black-blue-yellow-red-purple colormap
function col = bbyrp(handles)

%This value is taken from the vertical slider on song_gui; appears to range
%from 0 to 1. Dummy value used for now.
%val = 1-get(handles.slide_colormap,'value');
val = .3;

check_derive = 1;

handles.bgcolor = [.5 .5 .5];

%Dummy value used since, I think, this is always 1
%if get(handles.check_derivative,'value')==1
if check_derive == 1    
    col = repmat(linspace(0,1,256)',1,3);
    %set(handles.axes_spectrum,'color',handles.bgcolor);
    set(handles.axesSpecgram,'color',handles.bgcolor);
    coeffs = [.5 .5 .5]-handles.bgcolor;
    col(:,1) = col(:,1) - coeffs(1)*exp(-5*linspace(-1,1,256).^2)';
    col(:,2) = col(:,2) - coeffs(2)*exp(-5*linspace(-1,1,256).^2)';
    col(:,3) = col(:,3) - coeffs(3)*exp(-5*linspace(-1,1,256).^2)';
    col(find(col<0)) = 0;
    col(find(col>1)) = 1;
    set(handles.axesSpecgram,'clim',[-pi/2 pi/2]*(val+eps)^7);
else
    numblack = round(100*0.75*val);
    numother = round(100*(1-0.75*val)/4);

    col1 = zeros(numblack,3); % Black
    col2 = [linspace(0,0,numother)' linspace(0,0,numother)' linspace(0,0.5,numother)']; % Black to blue
    col3 = [linspace(0,1,numother)' linspace(0,1,numother)' linspace(0.5,0,numother)']; % Blue to yellow
    col4 = [linspace(1,1,numother)' linspace(1,0,numother)' linspace(0,0,numother)']; % Yellow to red
    col5 = [linspace(1,0.25,numother)' linspace(0,0,numother)' linspace(0,0.25,numother)']; % Red to purple

    col = [col1; col2; col3; col4; col5];
end



function editEdgeThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to editEdgeThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editEdgeThreshold as text
%        str2double(get(hObject,'String')) returns contents of editEdgeThreshold as a double
handles.edges=str2double(get(hObject,'String'));
set(handles.editEdgeThreshold,'String',num2str(handles.edges));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function editEdgeThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editEdgeThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







% --- Executes on button press in checkboxEdgefind.
function checkboxEdgefind_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxEdgefind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxEdgefind





function editPowerThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to editPowerThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPowerThreshold as text
%        str2double(get(hObject,'String')) returns contents of editPowerThreshold as a double
handles.threshold=str2double(get(hObject,'String'));
set(handles.editPowerThreshold,'String',num2str(handles.threshold));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editPowerThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPowerThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editSyllAbs_Callback(hObject, eventdata, handles)
% hObject    handle to editSyllAbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSyllAbs as text
%        str2double(get(hObject,'String')) returns contents of editSyllAbs as a double

handles.thresSyll=str2double(get(hObject,'String'));
set(handles.editSyllAbs,'String',num2str(handles.thresSyll));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editSyllAbs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSyllAbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function editEdgeAbs_Callback(hObject, eventdata, handles)
% hObject    handle to editEdgeAbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editEdgeAbs as text
%        str2double(get(hObject,'String')) returns contents of editEdgeAbs as a double
handles.thresEdge=str2double(get(hObject,'String'));
set(handles.editEdgeAbs,'String',num2str(handles.thresEdge));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editEdgeAbs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editEdgeAbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in buttonHelp.
function buttonHelp_Callback(hObject, eventdata, handles)
% hObject    handle to buttonHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

help_string=['comma(,) - zoom back' sprintf('\n') 'period(.) - zoom forward' sprintf('\n') 'p - play sound' sprintf('\n') 'space bar - jump to next syllable' sprintf('\n') sprintf('\n') 'Labeling Syllables:' sprintf('\n') sprintf('\n') '1-10 - labels syllables 1-10' sprintf('\n') 'u - unknown syllable' sprintf('\n') 'c - call' sprintf('\n') 's - subsong' sprintf('\n') 'delete - delete syllable' sprintf('\n') 'm - delete pause'];
handles.help=helpdlg(help_string,'Key to keys');

guidata(hObject, handles);



% --- Executes on button press in buttonLoadCell.
function buttonLoadCell_Callback(hObject, eventdata, handles)
% hObject    handle to buttonLoadCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.cellName,handles.cellPath]=uigetfile('*.mat');
load([handles.cellPath, handles.cellName]);
if (strcmp(presentCell(1).birdname,handles.exper.birdname) || presentCell(1).channel==handles.chan || presentCell(1).cell_no==handles.cell_no)
    if ~isfield(presentCell,'spont')
        for i=1:size(presentCell,2) %make the old cell files up to date with the spont field
            presentCell(i).spont=[];
        end
    end
   handles.cell_allrecs=presentCell;
   filenum = getFieldVector(presentCell,'filenum');
  
   index=find(filenum==handles.filenum);
   if isempty(index)
       
   else
       handles.cell=presentCell(index);
   end
    % see if there is data for the present filenu

   handles=ShowSpikes(handles);
else
    errordlg('loading a cell different from what you are working on')
    return;
end
guidata(hObject, handles);



% --- Executes on button press in buttonNewCell.
function buttonNewCell_Callback(hObject, eventdata, handles)
% hObject    handle to buttonNewCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%this will make a new cell. if a cell is laready there that cell will be
%abandoned at this point, and has to be reloaded to continue.
handles.cell.birdname=handles.exper.birdname;
handles.cell.filename=handles.spike_filename;
handles.cell.channel=handles.chan;
handles.cell.cell_no=handles.cell_no;
handles.cell.filenum=handles.filenum;
handles.cell.spikes=[];
handles.cell.noise=[];
handles.cell.spont=[];
[handles.cellName,handles.cellPath]=uiputfile('cell.mat','Save cell as:');
presentCell(1)=handles.cell; 
handles.cell_allrecs=presentCell;
save ([handles.cellPath,handles.cellName], 'presentCell');
guidata(hObject, handles);


% --- Executes on button press in buttonAddSpike.
function buttonAddSpike_Callback(hObject, eventdata, handles)
% hObject    handle to buttonAddSpike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (~isempty(handles.spiketime) && ~isempty(handles.unsorted_spikes))
    added_spike1 = find (handles.unsorted_spikes<handles.spiketime, 1,'last');
    added_spike2 = find (handles.unsorted_spikes>handles.spiketime, 1,'first');
    if isempty(added_spike2)
        added_spike=added_spike1;
    elseif (abs(handles.unsorted_spikes(added_spike1)-handles.spiketime)<abs(handles.unsorted_spikes(added_spike2)-handles.spiketime))
        added_spike=added_spike1;
    else
        added_spike=added_spike2;
    end
    added_spike_index=find(handles.sorted_spikes<handles.unsorted_spikes(added_spike), 1,'last');
    handles.sorted_spikes=[handles.sorted_spikes(1:added_spike_index); handles.unsorted_spikes(added_spike); handles.sorted_spikes((added_spike_index+1):end)];
    handles.unsorted_spikes=[handles.unsorted_spikes(1:added_spike-1); handles.unsorted_spikes(added_spike+1:end)];
    handles.spiketime=[];
    subplot(handles.axesRasters);
%     xlimit=get(gca,'xlim');
    cla;
% axis off;
     for i=1:length(handles.sorted_spikes)
         line([handles.sorted_spikes(i) handles.sorted_spikes(i)],[0.3 1],'color','r');
     end
    hold on;
    y = zeros(length(handles.unsorted_spikes));
    plot(handles.unsorted_spikes,y,'b+');
    hold off;
    handles.spikemarker=[];
    handles.cell.spikes=handles.sorted_spikes;
    handles.cell.noise=handles.unsorted_spikes;
end
guidata(hObject, handles);
% --- Executes on button press in buttonDeleteSpike.

function buttonDeleteSpike_Callback(hObject, eventdata, handles)
% hObject    handle to buttonDeleteSpike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (~isempty(handles.spiketime))
    added_spike1 = find (handles.sorted_spikes<handles.spiketime, 1,'last');
    added_spike2 = find (handles.sorted_spikes>handles.spiketime, 1,'first');
    if (abs(handles.sorted_spikes(added_spike1)-handles.spiketime)<abs(handles.sorted_spikes(added_spike2)-handles.spiketime));
        added_spike=added_spike1;
    else
        added_spike=added_spike2;
    end
    added_spike_index=find(handles.unsorted_spikes<handles.sorted_spikes(added_spike), 1,'last');
    if isempty(added_spike_index)
      added_spike_index=1;
    end
    if isempty(handles.unsorted_spikes)
       handles.unsorted_spikes=handles.sorted_spikes(added_spike);
    else
        handles.unsorted_spikes=[handles.unsorted_spikes(1:added_spike_index); handles.sorted_spikes(added_spike); handles.unsorted_spikes((added_spike_index+1):end)];
    end
    handles.sorted_spikes=[handles.sorted_spikes(1:added_spike-1); handles.sorted_spikes(added_spike+1:end)];
    handles.spiketime=[];
    
    subplot(handles.axesRasters);
%     xlimit=get(gca,'xlim');
    cla;
% axis off;
     for i=1:length(handles.sorted_spikes)
         line([handles.sorted_spikes(i) handles.sorted_spikes(i)],[0.3 1],'color','r');
     end
    hold on;
    y = zeros(length(handles.unsorted_spikes));
    plot(handles.unsorted_spikes,y,'b+');
    hold off;
    
end
handles.cell.spikes=handles.sorted_spikes;
handles.cell.noise=handles.unsorted_spikes;
handles.spikemarker=[];
guidata(hObject, handles);


% --- Executes on button press in checkboxShowSpikes.
function checkboxShowSpikes_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxShowSpikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxShowSpikes
handles=ShowSpikes(handles);
guidata(hObject, handles);

function handles=ShowSpikes(handles)

axes(handles.axesRasters);
cla;
if (get(handles.checkboxShowSpikes,'Value') == 1)
    % Checkbox is checked-take approriate action
    %see if there is a cell loaded
    if (isfield(handles,'cell_allrecs'))
       fieldVector=getFieldVector(handles.cell_allrecs,'filenum');
       
       index=find(fieldVector==handles.filenum);
        %and if it has the spikes for the filenumber and cell in question
        if (~isempty(index)&& handles.cell.channel==handles.chan && handles.cell.cell_no==handles.cell_no && ~isempty(handles.cell_allrecs(index).spikes))
            handles.sorted_spikes=handles.cell_allrecs(index).spikes;
            handles.unsorted_spikes=handles.cell_allrecs(index).noise;
            handles.cell.spikes=handles.sorted_spikes;
            handles.cell.noise=handles.unsorted_spikes;
             if (isempty(handles.cell.spikes))
                warndlg('no spikes!');
                uiwait;
                return;
            end
            axes(handles.axesRasters);
            for i=1:length(handles.sorted_spikes)
                line([handles.sorted_spikes(i) handles.sorted_spikes(i)],[0.3 1],'color','k');
            end
            hold on;
            y = zeros(length(handles.unsorted_spikes));
            plot(handles.unsorted_spikes,y,'b+');ylim([-0.3 1]);
            hold off;
        else 
            handles=loadSpikes(handles);
        end
    else
        handles=loadSpikes(handles); 
     end
else
    % Checkbox is not checked-take appropriate action
end

function handles=loadSpikes(handles)

handles.spike_filename = getspikeDatafile(handles.exper,handles.filenum,handles.chan);
handles.spikemarker=[];
if (handles.spike_filename=='0')
    warndlg('no spikes for this file');
    uiwait;
else
    load([handles.spike_filename]);
    handles.cell.filename=handles.spike_filename;
    cluster_class(:,2)=cluster_class(:,2)/1000; %threshhold crossings in seconds
    popup_string=(1:1:(max(cluster_class(:,1))));
    set(handles.popupCell,'String',popup_string);
    set(handles.popupCell,'Value',handles.cell_no);
    cellindex=find(cluster_class(:,1)==handles.cell_no);
    handles.sorted_spikes=cluster_class(cellindex,2);
    handles.cell.spikes=handles.sorted_spikes;
    cellindex=find(cluster_class(:,1)~=handles.cell_no); %fix this
    cellindex=sort(cellindex);
    handles.unsorted_spikes=cluster_class(cellindex,2);
    handles.cell.noise=handles.unsorted_spikes;
    handles.cell.spont=[];
    axes(handles.axesRasters);
    for i=1:length(handles.sorted_spikes)
        line([handles.sorted_spikes(i) handles.sorted_spikes(i)],[0.3 1],'color','r');
    end
    hold on;
    y = zeros(length(handles.unsorted_spikes));
    plot(handles.unsorted_spikes,y,'b+');ylim([-0.3 1]);
    hold off;
end

% --- Executes on button press in checkboxAutoID.
function checkboxAutoID_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAutoID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxAutoID



% --- Executes on selection change in popupCell.
function popupCell_Callback(hObject, eventdata, handles)
% hObject    handle to popupCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupCell contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupCell
handles.cell_no=get(hObject,'Value');
handles=ShowSpikes(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupCell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupDrug_Callback(hObject, eventdata, handles)
% hObject    handle to popupDrug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupDrug contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupDrug

contents = get(hObject,'String');% returns popupDrug contents as cell array
handles.drugstatus=contents{get(hObject,'Value')};
handles.drugindex=get(hObject,'Value');
if(handles.audioAnnotation.containsKey(handles.filename))
    currAnnot = handles.audioAnnotation.get(handles.filename);
    currAnnot.drugindex=handles.drugindex;
    currAnnot.drugstatus=handles.drugstatus;
    handles.audioAnnotation.put(handles.filename, currAnnot);
end

guidata(hObject, handles);

function popupDrug_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupDrug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function buttonSpont_Callback(hObject, eventdata, handles)
% hObject    handle to buttonSpont (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (isempty(handles.cell.spikes))
    warndlg('no spikes!');
    uiwait;
    return;
else
    startTime = (handles.startNdx-1)/handles.fs;
    endTime = (handles.endNdx-1)/handles.fs;
    indexSpikes=find(handles.cell.spikes>startTime & handles.cell.spikes<endTime);
    spontSpikes=handles.cell.spikes(indexSpikes);
    handles.cell.spont=diff(spontSpikes);
    
end
guidata(hObject, handles);

function push_testsono_Callback(hObject, eventdata, handles)
% hObject    handle to push_testsono (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axesSpecgram);
cla
hold on

colormap(bbyrp(handles))
handles.spectrumplots = [];
handles.wav = handles.audio;
handles.bgcolor = [.5 .5 .5];
%This whole top section is not normally used in song_gui and I don't think
%relevant here
% if get(handles.check_usesegments,'value')==1
%     seg_list = find(handles.selectedseg==1);
%     if get(handles.check_currentwindow,'value')==1
%         isin = (handles.segments/handles.fs > handles.starttime) & (handles.segments/handles.fs < handles.starttime+str2num(get(handles.edit_timescale,'string')));
%         isin = sum(isin,2)';
%         seg_list = intersect(seg_list,find(isin>0));
%     end
%     for c = seg_list
%         [p freq t] = power_spectrum(handles.wav(max([1 round(handles.segments(c,1))]):min([round(handles.segments(c,2)) length(handles.wav)])),handles.fs,get(handles.check_derivative,'value'));
%         handles.satur = str2num(get(handles.edit_saturation,'string'));
%         if get(handles.check_derivative,'value')==1
%             set(gca,'color',handles.bgcolor);
%             subplot(handles.axes_spectrum);
%             h = imagesc(t+handles.segments(c,1)/handles.fs,freq/1000,atan(p/(handles.satur*10^10)));
%         else
%             set(gca,'color',[0 0 0]);
%             subplot(handles.axes_spectrum);
%             h = imagesc(t+handles.segments(c,1)/handles.fs,freq/1000,p);
%             set(gca,'clim',[0 handles.satur]);
%         end
%         set(h,'ButtonDownFcn','song_gui(''clicksonogram'',gcbo,[],guidata(gcbo))');
%         drawnow
%     end
%else
%     if get(handles.check_currentwindow,'value')==0
%         lst = [0:handles.fs:length(handles.wav) length(handles.wav)];
%    else
%          %lst =fix(handles.starttime)*handles.fs:handles.fs:ceil(handles.starttime+str2num(get(handles.edit_timescale,'string')))*handles.fs;
%          lst = fix(handles.startNdx):handles.fs:ceil(handles.endNdx);
%          lst = lst(find(lst<length(handles.wav)));
%          if lst(end) > length(handles.wav)-handles.fs
%              lst = [lst length(handles.wav)];
%          end
%    end
%    for c = 1:length(lst)-1
        %[p freq t] = power_spectrum(handles.wav(lst(c)+1:min([lst(c+1)+round(0.1*handles.fs) length(handles.wav)])),handles.fs,get(handles.check_derivative,'value'));
        %handles.satur = str2num(get(handles.edit_saturation,'string'));
        %[p freq t] = power_spectrum(handles.wav(lst(c)+1:min([lst(c+1)+round(0.1*handles.fs) length(handles.wav)])),handles.fs,1);
        [p freq t] = power_spectrum(handles.wav(handles.startNdx:handles.endNdx),handles.fs,1);
        %[p freq t] = power_spectrum(handles.wav,handles.fs,1);
        handles.satur = 25; %looks like this number was always 25
        %if get(handles.check_derivative,'value')==1
            set(gca,'color',handles.bgcolor);
            axes(handles.axesSpecgram);
            %h = imagesc(t+lst(c)/handles.fs,freq/1000,atan(p/(handles.satur*10^10)));
            h = imagesc(t,freq/1000,atan(p/(handles.satur*10^10)));
%         else
%             set(gca,'color',[0 0 0]);
%             subplot(handles.axes_spectrum);
%             h = imagesc(t+lst(c)/handles.fs,freq/1000,p);
%             set(gca,'clim',[0 handles.satur]);
%         end
        %set(h,'ButtonDownFcn','song_gui(''cb_specgram_click'',gcbo,[],guidata(gcbo))');
        set(h, 'ButtonDownFcn', @cb_specgram_click);
        drawnow
%    end
%end
    
%subplot(handles.axes_spectrum);
subplot(handles.axesSpecgram);
set(gca,'ydir','normal');
ylim([.5 9]);
hold off
colormap(bbyrp(handles))
xlabel('Time (sec)');
ylabel('Frequency (kHz)')
%set(handles.text_warn,'visible','off')


% --- Power spectrum of a song segment
function [p f t] = power_spectrum(wv,fs,issd)

winstep = 0.1;
npad = 1024;

E = [0.10082636028528   0.48810350894928
   0.10532280057669   0.50016260147095
   0.10990284383297   0.51226782798767
   0.11456660926342   0.52441596984863
   0.11931420862675   0.53660380840302
   0.12414572387934   0.54882800579071
   0.12906123697758   0.56108516454697
   0.13406077027321   0.57337194681168
   0.13914434611797   0.58568495512009
   0.14431197941303   0.59802067279816
   0.14956362545490   0.61037564277649
   0.15489923954010   0.62274640798569
   0.16031874716282   0.63512933254242
   0.16582205891609   0.64752095937729
   0.17140907049179   0.65991759300232
   0.17707960307598   0.67231565713882
   0.18283350765705   0.68471157550812
   0.18867060542107   0.69710153341293
   0.19459065794945   0.70948195457458
   0.20059344172478   0.72184902429581
   0.20667870342732   0.73419916629791
   0.21284613013268   0.74652844667435
   0.21909542381764   0.75883322954178
   0.22542624175549   0.77110970020294
   0.23183822631836   0.78335398435593
   0.23833100497723   0.79556238651276
   0.24490414559841   0.80773097276688
   0.25155723094940   0.81985598802566
   0.25828978419304   0.83193349838257
   0.26510134339333   0.84395974874497
   0.27199140191078   0.85593080520630
   0.27895939350128   0.86784279346466
   0.28600478172302   0.87969189882278
   0.29312700033188   0.89147418737411
   0.30032539367676   0.90318578481674
   0.30759936571121   0.91482281684875
   0.31494826078415   0.92638140916824
   0.32237139344215   0.93785762786865
   0.32986804842949   0.94924771785736
   0.33743748068809   0.96054774522781
   0.34507897496223   0.97175377607346
   0.35279172658920   0.98286205530167
   0.36057496070862   0.99386870861054
   0.36842781305313   1.00476980209351
   0.37634941935539   1.01556169986725
   0.38433897495270   1.02624034881592
   0.39239549636841   1.03680217266083
   0.40051811933518   1.04724323749542
   0.40870589017868   1.05755984783173
   0.41695779561996   1.06774818897247
   0.42527288198471   1.07780468463898
   0.43365013599396   1.08772540092468
   0.44208851456642   1.09750676155090
   0.45058691501617   1.10714507102966
   0.45914432406425   1.11663687229156
   0.46775957942009   1.12597823143005
   0.47643154859543   1.13516581058502
   0.48515912890434   1.14419603347778
   0.49394109845161   1.15306520462036
   0.50277632474899   1.16176998615265
   0.51166349649429   1.17030680179596
   0.52060145139694   1.17867243289948
   0.52958887815475   1.18686318397522
   0.53862458467483   1.19487595558167
   0.54770720005035   1.20270740985870
   0.55683541297913   1.21035408973694
   0.56600785255432   1.21781289577484
   0.57522326707840   1.22508060932159
   0.58448016643524   1.23215413093567
   0.59377723932266   1.23903024196625
   0.60311305522919   1.24570596218109
   0.61248612403870   1.25217819213867
   0.62189501523972   1.25844407081604
   0.63133829832077   1.26450049877167
   0.64081448316574   1.27034485340118
   0.65032202005386   1.27597427368164
   0.65985941886902   1.28138577938080
   0.66942512989044   1.28657686710358
   0.67901766300201   1.29154479503632
   0.68863534927368   1.29628694057465
   0.69827669858932   1.30080080032349
   0.70794004201889   1.30508399009705
   0.71762382984161   1.30913388729095
   0.72732639312744   1.31294822692871
   0.73704612255096   1.31652474403381
   0.74678134918213   1.31986105442047
   0.75653040409088   1.32295513153076
   0.76629155874252   1.32580471038818
   0.77606326341629   1.32840788364410
   0.78584367036820   1.33076262474060
   0.79563117027283   1.33286690711975
   0.80542397499084   1.33471894264221
   0.81522035598755   1.33631694316864
   0.82501864433289   1.33765923976898
   0.83481699228287   1.33874404430389
   0.84461367130279   1.33956992626190
   0.85440695285797   1.34013533592224
   0.86419498920441   1.34043884277344
   0.87397605180740   1.34047901630402
   0.88374835252762   1.34025466442108
   0.89351004362106   1.33976447582245
   0.90325933694840   1.33900737762451
   0.91299438476563   1.33798217773438
   0.92271345853805   1.33668816089630
   0.93241471052170   1.33512413501740
   0.94209629297257   1.33328926563263
   0.95175635814667   1.33118307590485
   0.96139311790466   1.32880449295044
   0.97100472450256   1.32615327835083
   0.98058933019638   1.32322859764099
   0.99014514684677   1.32003021240234
   0.99967026710510   1.31655764579773
   1.00916278362274   1.31281065940857
   1.01862108707428   1.30878901481628
   1.02804315090179   1.30449247360229
   1.03742706775665   1.29992127418518
   1.04677128791809   1.29507517814636
   1.05607366561890   1.28995430469513
   1.06533253192902   1.28455901145935
   1.07454609870911   1.27888941764832
   1.08371245861053   1.27294600009918
   1.09282970428467   1.26672899723053
   1.10189616680145   1.26023912429810
   1.11090993881226   1.25347697734833
   1.11986923217773   1.24644303321838
   1.12877237796783   1.23913824558258
   1.13761734962463   1.23156344890594
   1.14640247821808   1.22371935844421
   1.15512585639954   1.21560716629028
   1.16378593444824   1.20722794532776
   1.17238080501556   1.19858264923096
   1.18090879917145   1.18967282772064
   1.18936800956726   1.18049955368042
   1.19775688648224   1.17106437683105
   1.20607352256775   1.16136872768402
   1.21431636810303   1.15141415596008
   1.22248351573944   1.14120233058929
   1.23057353496552   1.13073492050171
   1.23858451843262   1.12001371383667
   1.24651479721069   1.10904061794281
   1.25436294078827   1.09781765937805
   1.26212716102600   1.08634674549103
   1.26980578899384   1.07462990283966
   1.27739727497101   1.06266951560974
   1.28489995002747   1.05046772956848
   1.29231238365173   1.03802680969238
   1.29963290691376   1.02534937858582
   1.30685997009277   1.01243758201599
   1.31399202346802   0.99929422140122
   1.32102763652802   0.98592185974121
   1.32796525955200   0.97232311964035
   1.33480334281921   0.95850080251694
   1.34154057502747   0.94445770978928
   1.34817540645599   0.93019676208496
   1.35470652580261   0.91572093963623
   1.36113238334656   0.90103322267532
   1.36745178699493   0.88613671064377
   1.37366318702698   0.87103462219238
   1.37976539134979   0.85573017597198
   1.38575696945190   0.84022665023804
   1.39163672924042   0.82452732324600
   1.39740324020386   0.80863571166992
   1.40305554866791   0.79255521297455
   1.40859210491180   0.77628940343857
   1.41401195526123   0.75984185934067
   1.41931378841400   0.74321621656418
   1.42449653148651   0.72641617059708
   1.42955899238586   0.70944553613663
   1.43450009822845   0.69230806827545
   1.43931877613068   0.67500764131546
   1.44401395320892   0.65754824876785
   1.44858467578888   0.63993370532990
   1.45302975177765   0.62216812372208
   1.45734846591949   0.60425555706024
   1.46153962612152   0.58620011806488
   1.46560251712799   0.56800597906113
   1.46953618526459   0.54967725276947
   1.47333967685699   0.53121829032898
   1.47701227664948   0.51263326406479
   1.48055315017700   0.49392655491829
   1.48396146297455   0.47510254383087
   1.48723638057709   0.45616558194161
   1.49037742614746   0.43712010979652
   1.49338376522064   0.41797062754631
   1.49625468254089   0.39872157573700
   1.49898970127106   0.37937757372856
   1.50158810615540   0.35994309186935
   1.50404930114746   0.34042277932167
   1.50637280941010   0.32082122564316
   1.50855803489685   0.30114310979843
   1.51060461997986   0.28139305114746
   1.51251196861267   0.26157575845718
   1.51427984237671   0.24169597029686
   1.51590764522552   0.22175838053226
   1.51739513874054   0.20176777243614
   1.51874196529388   0.18172888457775
   1.51994776725769   0.16164653003216
   1.52101242542267   0.14152547717094
   1.52193558216095   0.12137054651976
   1.52271711826324   0.10118655860424
   1.52335667610168   0.08097834140062
   1.52385437488556   0.06075073406100
   1.52420985698700   0.04050857573748
   1.52442324161530   0.02025671303272
   1.52449429035187  -0.00000000211473
   1.52442324161530  -0.02025671675801
   1.52420985698700  -0.04050857946277
   1.52385437488556  -0.06075073778629
   1.52335667610168  -0.08097834885120
   1.52271711826324  -0.10118656605482
   1.52193558216095  -0.12137055397034
   1.52101242542267  -0.14152547717094
   1.51994776725769  -0.16164653003216
   1.51874196529388  -0.18172889947891
   1.51739513874054  -0.20176777243614
   1.51590764522552  -0.22175839543343
   1.51427984237671  -0.24169597029686
   1.51251196861267  -0.26157575845718
   1.51060461997986  -0.28139305114746
   1.50855803489685  -0.30114310979843
   1.50637280941010  -0.32082122564316
   1.50404930114746  -0.34042277932167
   1.50158810615540  -0.35994309186935
   1.49898970127106  -0.37937757372856
   1.49625468254089  -0.39872160553932
   1.49338376522064  -0.41797062754631
   1.49037742614746  -0.43712010979652
   1.48723638057709  -0.45616558194161
   1.48396146297455  -0.47510254383087
   1.48055315017700  -0.49392655491829
   1.47701227664948  -0.51263326406479
   1.47333967685699  -0.53121829032898
   1.46953618526459  -0.54967725276947
   1.46560251712799  -0.56800597906113
   1.46153962612152  -0.58620011806488
   1.45734846591949  -0.60425561666489
   1.45302975177765  -0.62216812372208
   1.44858467578888  -0.63993370532990
   1.44401395320892  -0.65754824876785
   1.43931877613068  -0.67500770092010
   1.43450009822845  -0.69230806827545
   1.42955899238586  -0.70944553613663
   1.42449653148651  -0.72641617059708
   1.41931378841400  -0.74321621656418
   1.41401195526123  -0.75984185934067
   1.40859210491180  -0.77628940343857
   1.40305554866791  -0.79255521297455
   1.39740324020386  -0.80863571166992
   1.39163672924042  -0.82452732324600
   1.38575696945190  -0.84022665023804
   1.37976539134979  -0.85573017597198
   1.37366318702698  -0.87103468179703
   1.36745178699493  -0.88613677024841
   1.36113238334656  -0.90103322267532
   1.35470652580261  -0.91572093963623
   1.34817540645599  -0.93019676208496
   1.34154057502747  -0.94445770978928
   1.33480334281921  -0.95850080251694
   1.32796525955200  -0.97232311964035
   1.32102763652802  -0.98592185974121
   1.31399202346802  -0.99929428100586
   1.30685997009277  -1.01243758201599
   1.29963290691376  -1.02534937858582
   1.29231238365173  -1.03802680969238
   1.28489995002747  -1.05046772956848
   1.27739727497101  -1.06266951560974
   1.26980578899384  -1.07462990283966
   1.26212716102600  -1.08634674549103
   1.25436294078827  -1.09781765937805
   1.24651479721069  -1.10904061794281
   1.23858451843262  -1.12001371383667
   1.23057353496552  -1.13073492050171
   1.22248351573944  -1.14120233058929
   1.21431636810303  -1.15141415596008
   1.20607352256775  -1.16136872768402
   1.19775688648224  -1.17106437683105
   1.18936800956726  -1.18049955368042
   1.18090879917145  -1.18967282772064
   1.17238080501556  -1.19858264923096
   1.16378593444824  -1.20722794532776
   1.15512585639954  -1.21560716629028
   1.14640247821808  -1.22371935844421
   1.13761734962463  -1.23156344890594
   1.12877237796783  -1.23913824558258
   1.11986923217773  -1.24644303321838
   1.11090993881226  -1.25347697734833
   1.10189616680145  -1.26023912429810
   1.09282970428467  -1.26672899723053
   1.08371245861053  -1.27294600009918
   1.07454609870911  -1.27888941764832
   1.06533253192902  -1.28455901145935
   1.05607366561890  -1.28995430469513
   1.04677128791809  -1.29507517814636
   1.03742706775665  -1.29992127418518
   1.02804315090179  -1.30449247360229
   1.01862108707428  -1.30878901481628
   1.00916278362274  -1.31281065940857
   0.99967020750046  -1.31655764579773
   0.99014514684677  -1.32003021240234
   0.98058933019638  -1.32322859764099
   0.97100472450256  -1.32615327835083
   0.96139311790466  -1.32880449295044
   0.95175635814667  -1.33118307590485
   0.94209629297257  -1.33328926563263
   0.93241471052170  -1.33512413501740
   0.92271345853805  -1.33668816089630
   0.91299438476563  -1.33798217773438
   0.90325933694840  -1.33900737762451
   0.89351004362106  -1.33976447582245
   0.88374835252762  -1.34025466442108
   0.87397605180740  -1.34047901630402
   0.86419498920441  -1.34043884277344
   0.85440695285797  -1.34013533592224
   0.84461367130279  -1.33956992626190
   0.83481699228287  -1.33874404430389
   0.82501864433289  -1.33765923976898
   0.81522035598755  -1.33631694316864
   0.80542397499084  -1.33471894264221
   0.79563117027283  -1.33286690711975
   0.78584367036820  -1.33076262474060
   0.77606326341629  -1.32840788364410
   0.76629155874252  -1.32580471038818
   0.75653040409088  -1.32295513153076
   0.74678134918213  -1.31986105442047
   0.73704612255096  -1.31652474403381
   0.72732639312744  -1.31294822692871
   0.71762382984161  -1.30913388729095
   0.70794004201889  -1.30508399009705
   0.69827669858932  -1.30080080032349
   0.68863534927368  -1.29628694057465
   0.67901766300201  -1.29154479503632
   0.66942512989044  -1.28657686710358
   0.65985941886902  -1.28138577938080
   0.65032202005386  -1.27597427368164
   0.64081448316574  -1.27034485340118
   0.63133829832077  -1.26450049877167
   0.62189501523972  -1.25844407081604
   0.61248612403870  -1.25217819213867
   0.60311305522919  -1.24570596218109
   0.59377723932266  -1.23903024196625
   0.58448016643524  -1.23215413093567
   0.57522326707840  -1.22508060932159
   0.56600785255432  -1.21781289577484
   0.55683541297913  -1.21035408973694
   0.54770720005035  -1.20270740985870
   0.53862458467483  -1.19487595558167
   0.52958887815475  -1.18686318397522
   0.52060145139694  -1.17867243289948
   0.51166349649429  -1.17030680179596
   0.50277632474899  -1.16176998615265
   0.49394109845161  -1.15306520462036
   0.48515912890434  -1.14419603347778
   0.47643154859543  -1.13516581058502
   0.46775957942009  -1.12597823143005
   0.45914432406425  -1.11663687229156
   0.45058691501617  -1.10714507102966
   0.44208851456642  -1.09750676155090
   0.43365013599396  -1.08772540092468
   0.42527288198471  -1.07780468463898
   0.41695779561996  -1.06774818897247
   0.40870586037636  -1.05755984783173
   0.40051811933518  -1.04724323749542
   0.39239549636841  -1.03680217266083
   0.38433894515037  -1.02624034881592
   0.37634941935539  -1.01556169986725
   0.36842781305313  -1.00476980209351
   0.36057496070862  -0.99386870861054
   0.35279172658920  -0.98286205530167
   0.34507897496223  -0.97175377607346
   0.33743748068809  -0.96054774522781
   0.32986804842949  -0.94924771785736
   0.32237139344215  -0.93785762786865
   0.31494826078415  -0.92638140916824
   0.30759936571121  -0.91482281684875
   0.30032539367676  -0.90318578481674
   0.29312697052956  -0.89147418737411
   0.28600478172302  -0.87969189882278
   0.27895939350128  -0.86784279346466
   0.27199140191078  -0.85593080520630
   0.26510134339333  -0.84395974874497
   0.25828978419304  -0.83193349838257
   0.25155723094940  -0.81985598802566
   0.24490414559841  -0.80773097276688
   0.23833100497723  -0.79556238651276
   0.23183822631836  -0.78335398435593
   0.22542624175549  -0.77110970020294
   0.21909542381764  -0.75883322954178
   0.21284613013268  -0.74652844667435
   0.20667870342732  -0.73419916629791
   0.20059344172478  -0.72184902429581
   0.19459065794945  -0.70948195457458
   0.18867060542107  -0.69710153341293
   0.18283350765705  -0.68471157550812
   0.17707960307598  -0.67231565713882
   0.17140905559063  -0.65991759300232
   0.16582205891609  -0.64752095937729
   0.16031874716282  -0.63512933254242
   0.15489923954010  -0.62274640798569
   0.14956361055374  -0.61037564277649
   0.14431196451187  -0.59802067279816
   0.13914434611797  -0.58568495512009
   0.13406075537205  -0.57337194681168
   0.12906122207642  -0.56108516454697
   0.12414572387934  -0.54882800579071
   0.11931420117617  -0.53660380840302
   0.11456660181284  -0.52441596984863
   0.10990284383297  -0.51226782798767
   0.10532280057669  -0.50016260147095
   0.10082635283470  -0.48810350894928];

winlen = size(E,1);
winstep = round(winlen*winstep);

[x y] = meshgrid(1:winstep:length(wv)-winlen,0:winlen-1);
vals = wv(x+y);
clear x y
TSM=vals';

multconst = 27539;
J1=(fft(TSM(:,:).*(ones(size(TSM,1),1)*(E(:,1))'),2^10,2)) * multconst;
J2=(fft(TSM(:,:).*(ones(size(TSM,1),1)*(E(:,2))'),2^10,2)) * multconst;

if issd==1
    m_time_deriv=-1*(real(J1).*real(J2)+imag(J1).*imag(J2));
    m_freq_deriv=((imag(J1).*real(J2)-real(J1).*imag(J2)));
    m_time_deriv_max=max(m_time_deriv.^2,[],2);
    m_freq_deriv_max=max(m_freq_deriv.^2,[],2);
    m_FM=atan(m_time_deriv_max./(m_freq_deriv_max+eps));
    cFM=cos(m_FM);
    sFM=sin(m_FM);
    m_powSpec = m_time_deriv.*(sFM*ones(1,size(m_time_deriv,2)))+m_freq_deriv.*(cFM*ones(1,size(m_freq_deriv,2)));
else
    m_powSpec = real(J1).^2+real(J2).^2+imag(J1).^2+imag(J1).^2;
end

if issd==1
    p = m_powSpec;
else
    p = log(m_powSpec+eps);
end

t = (0:size(m_powSpec,1)-1)*winstep/fs;
f = fs*(0:size(m_powSpec,2)-1)/size(m_powSpec,2);

num = size(p,2)/2;
p = p(:,1:num)';
f = f(1:num);
