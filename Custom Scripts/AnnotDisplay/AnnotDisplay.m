function varargout = AnnotDisplay(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AnnotDisplay_OpeningFcn, ...
                   'gui_OutputFcn',  @AnnotDisplay_OutputFcn, ...
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
function AnnotDisplay_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
function varargout = AnnotDisplay_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in push_loadannot.
function push_loadannot_Callback(hObject, eventdata, handles)

[file,path] = uigetfile('*.mat', 'Create a .mat file for the audioAnnotation:');

if ~(isequal(file,0) || isequal(path,0))
    
    % LOAD ANNOTFILE
    handles.annotFilename = [path file];
    load(handles.annotFilename);
    sylStart=[];
    sylEnd=[];
    sylType=[]; 
    IndxKey=[];
    
    for i=1:length(elements)
        indx_i=repmat(i,length(elements{i}.segType),1);
        IndxKey=vertcat(IndxKey,indx_i);
        sylStart=horzcat(sylStart, elements{i}.segFileStartTimes);
        sylEnd=horzcat(sylEnd, elements{i}.segFileEndTimes);
        sylType=vertcat(sylType, elements{i}.segType);
    end
    
    % LOAD AUDIO
    audio=cell(length(keys),1);
    for i=1:length(keys)
        %audio{i} = audioread(keys{i})
        
        n = char(keys{i});
        % Strip the channels out of the raw recording file
        if strcmp(n(end-3:end), '.wav')
            %[r, fs] = wavread(keys{i});
            [r, fs] = audioread(keys{i});
            if size(r,2)>1
                audio{i} = r(:,1); %Take only the first channel from stereo audio recordings
            else
                audio{i} = r;
            end
        elseif strcmp(n(end-3:end), '.dat')
            [r, fs] = getChannels(keys{i});
            audio{i} = r(1,:)'; %Take only the first channel from multiplexed A&N recordings
        end
    end
    
    % PASTE IN HANDLES
    handles.audio=audio;
    handles.elements = elements;
    handles.keys = keys;
    handles.sylStart = sylStart*44150;
    handles.sylEnd = sylEnd*44150;
    handles.sylType = sylType;
    handles.IndxKey=IndxKey;
end

% PARAMETER FOR THE DISPLAY
set(handles.popupmenu,'Value',1)
handles.sylvalue=1;
handles.syltempIndx=find(handles.sylType==handles.sylvalue);
handles.syltempBin=handles.sylType==handles.sylvalue;
handles.syltempStart=handles.sylStart(handles.syltempBin);
handles.syltempEnd=handles.sylEnd(handles.syltempBin);
handles.syltempType=handles.sylType(handles.syltempBin);
handles.syltempIndxKey=handles.IndxKey(handles.syltempBin);
handles.temp_max_page=floor(length(handles.syltempType)/84);

% PAGE NUMBER
handles.page=1;
set(handles.popupmenu_page,'String',num2str(((1:handles.temp_max_page+1))'));
set(handles.popupmenu_page,'Value',handles.page);

handles.axes=[handles.axes1 handles.axes2 handles.axes3 handles.axes4 handles.axes5 handles.axes6 handles.axes7 handles.axes8 handles.axes9 handles.axes10 handles.axes11 handles.axes12 handles.axes13 handles.axes14 handles.axes15 handles.axes16 handles.axes17 handles.axes18 handles.axes19 handles.axes20 handles.axes21 handles.axes22 handles.axes23 handles.axes24 handles.axes25 handles.axes26 handles.axes27 handles.axes28 handles.axes29 handles.axes30 handles.axes31 handles.axes32 handles.axes33 handles.axes34 handles.axes35 handles.axes36 handles.axes37 handles.axes38 handles.axes39 handles.axes40 handles.axes41 handles.axes42 handles.axes43 handles.axes44 handles.axes45 handles.axes46 handles.axes47 handles.axes48 handles.axes49 handles.axes50 handles.axes51 handles.axes52 handles.axes53 handles.axes54 handles.axes55 handles.axes56 handles.axes57 handles.axes58 handles.axes59 handles.axes60 handles.axes61 handles.axes62 handles.axes63 handles.axes64 handles.axes65 handles.axes66 handles.axes67 handles.axes68 handles.axes69 handles.axes70 handles.axes71 handles.axes72 handles.axes73 handles.axes74 handles.axes75 handles.axes76 handles.axes77 handles.axes78 handles.axes79 handles.axes80 handles.axes81 handles.axes82 handles.axes83 handles.axes84];
handles.edit=[handles.edit1 handles.edit2 handles.edit3 handles.edit4 handles.edit5 handles.edit6 handles.edit7 handles.edit8 handles.edit9 handles.edit10 handles.edit11 handles.edit12 handles.edit13 handles.edit14 handles.edit15 handles.edit16 handles.edit17 handles.edit18 handles.edit19 handles.edit20 handles.edit21 handles.edit22 handles.edit23 handles.edit24 handles.edit25 handles.edit26 handles.edit27 handles.edit28 handles.edit29 handles.edit30 handles.edit31 handles.edit32 handles.edit33 handles.edit34 handles.edit35 handles.edit36 handles.edit37 handles.edit38 handles.edit39 handles.edit40 handles.edit41 handles.edit42 handles.edit43 handles.edit44 handles.edit45 handles.edit46 handles.edit47 handles.edit48 handles.edit49 handles.edit50 handles.edit51 handles.edit52 handles.edit53 handles.edit54 handles.edit55 handles.edit56 handles.edit57 handles.edit58 handles.edit59 handles.edit60 handles.edit61 handles.edit62 handles.edit63 handles.edit64 handles.edit65 handles.edit66 handles.edit67 handles.edit68 handles.edit69 handles.edit70 handles.edit71 handles.edit72 handles.edit73 handles.edit74 handles.edit75 handles.edit76 handles.edit77 handles.edit78 handles.edit79 handles.edit80 handles.edit81 handles.edit82 handles.edit83 handles.edit84];


% SET TIME WINDOW in ms, except 0 which mean no margin
set(handles.slider1,'Value',200)
set(handles.slider1,'Min',0)
set(handles.slider1,'Max',300)

% DISPLAY THE AXES
for i=1:84
    if (handles.page-1)*84+i < length(handles.syltempType)
        set(handles.edit(i),'String',handles.syltempType((handles.page-1)*84+i))
        set(handles.edit(i), 'ButtonDownFcn',@cb_remove_file)
    else
        set(handles.edit(i),'String',0)
        set(handles.edit(i),'BackgroundColor',[1 1 1])
    end
    handles = updateAudAxes(handles,(handles.page-1)*84+i,handles.axes(i));
end

guidata(hObject, handles);

function cb_remove_file(hObject, evnt)
handles = guidata(hObject);
set(hObject,'String',str2double(get(handles.edit_rightClick,'String')))
set(hObject,'BackgroundColor',[1 0.6 0.6])
guidata(hObject, handles);

% --- Executes on button press in push_savetemplate.
function push_savetemplate_Callback(hObject, eventdata, handles)
if ~isfield(handles,'annotFilename')
    warndlg('Nice try, load an annotation file before');
else
    % save
    for i=1:84
        if str2double(get(handles.edit(i),'String')) ~= 0
            handles.syltempType((handles.page-1)*84+i)= str2double(get(handles.edit(i),'String'));
        end
    end
    handles.sylType(handles.syltempIndx)=handles.syltempType(1:length(handles.sylType(handles.syltempIndx)));
    
    %EXPORT
    
    %load annot
    load(handles.annotFilename);
    
    %remove 111
    delIndx=find(handles.sylType==111);
    handles.sylType(delIndx) = [];
    handles.sylStart(delIndx) = [];
    handles.sylEnd(delIndx) = [];
    handles.IndxKey(delIndx) = [];
    
    %save in annot
    uniq=unique(handles.IndxKey);
    newelements=cell(1,length(uniq));
    newkeys=cell(1,length(uniq));
    audio=cell(1,length(uniq));
    parfor i=1:length(uniq)
        j=uniq(i);
        newelements{i}.exper=elements{j}.exper;
        newelements{i}.filenum=elements{j}.filenum;
        newelements{i}.segAbsStartTimes= getFileTime(keys{j}) + (handles.sylStart(handles.IndxKey==j)/(44150*(24*60*60)));
        newelements{i}.segFileStartTimes=handles.sylStart(handles.IndxKey==j)/44150;
        newelements{i}.segFileEndTimes=handles.sylEnd(handles.IndxKey==j)/44150;
        newelements{i}.segType=handles.sylType(handles.IndxKey==j);
        newelements{i}.fs=elements{j}.fs;
        newelements{i}.drugstatus=elements{j}.drugstatus;
        newelements{i}.directstatus=elements{j}.directstatus;
        newkeys{i}=keys{j}
        audio{i}=handles.audio{j}
    end
    elements=newelements;
    keys=newkeys;
    save(handles.annotFilename,'elements','keys');
    handles.audio=audio;
    result_save='done'
    
    
    %REFRESH PAGE
    
    sylStart=[];
    sylEnd=[];
    sylType=[];
    IndxKey=[];
    
    for i=1:length(elements)
        indx_i=repmat(i,length(elements{i}.segType),1);
        IndxKey=vertcat(IndxKey,indx_i);
        sylStart=horzcat(sylStart, elements{i}.segFileStartTimes);
        sylEnd=horzcat(sylEnd, elements{i}.segFileEndTimes);
        sylType=vertcat(sylType, elements{i}.segType);
    end
    
    
    handles.elements = elements;
    handles.keys = keys;
    handles.sylStart = sylStart*44150;
    handles.sylEnd = sylEnd*44150;
    handles.sylType = sylType;
    handles.IndxKey=IndxKey;
    
    % PARAMETER FOR THE DISPLAY
    set(handles.popupmenu,'Value',1)
    handles.sylvalue=1;
    handles.syltempIndx=find(handles.sylType==handles.sylvalue);
    handles.syltempBin=handles.sylType==handles.sylvalue;
    handles.syltempStart=handles.sylStart(handles.syltempBin);
    handles.syltempEnd=handles.sylEnd(handles.syltempBin);
    handles.syltempType=handles.sylType(handles.syltempBin);
    handles.syltempIndxKey=handles.IndxKey(handles.syltempBin);
    handles.temp_max_page=floor(length(handles.syltempType)/84);
    
    % PAGE NUMBER
    handles.page=1;
    set(handles.popupmenu_page,'String',num2str(((1:handles.temp_max_page+1))'));
    set(handles.popupmenu_page,'Value',handles.page);
    
    handles.axes=[handles.axes1 handles.axes2 handles.axes3 handles.axes4 handles.axes5 handles.axes6 handles.axes7 handles.axes8 handles.axes9 handles.axes10 handles.axes11 handles.axes12 handles.axes13 handles.axes14 handles.axes15 handles.axes16 handles.axes17 handles.axes18 handles.axes19 handles.axes20 handles.axes21 handles.axes22 handles.axes23 handles.axes24 handles.axes25 handles.axes26 handles.axes27 handles.axes28 handles.axes29 handles.axes30 handles.axes31 handles.axes32 handles.axes33 handles.axes34 handles.axes35 handles.axes36 handles.axes37 handles.axes38 handles.axes39 handles.axes40 handles.axes41 handles.axes42 handles.axes43 handles.axes44 handles.axes45 handles.axes46 handles.axes47 handles.axes48 handles.axes49 handles.axes50 handles.axes51 handles.axes52 handles.axes53 handles.axes54 handles.axes55 handles.axes56 handles.axes57 handles.axes58 handles.axes59 handles.axes60 handles.axes61 handles.axes62 handles.axes63 handles.axes64 handles.axes65 handles.axes66 handles.axes67 handles.axes68 handles.axes69 handles.axes70 handles.axes71 handles.axes72 handles.axes73 handles.axes74 handles.axes75 handles.axes76 handles.axes77 handles.axes78 handles.axes79 handles.axes80 handles.axes81 handles.axes82 handles.axes83 handles.axes84];
    handles.edit=[handles.edit1 handles.edit2 handles.edit3 handles.edit4 handles.edit5 handles.edit6 handles.edit7 handles.edit8 handles.edit9 handles.edit10 handles.edit11 handles.edit12 handles.edit13 handles.edit14 handles.edit15 handles.edit16 handles.edit17 handles.edit18 handles.edit19 handles.edit20 handles.edit21 handles.edit22 handles.edit23 handles.edit24 handles.edit25 handles.edit26 handles.edit27 handles.edit28 handles.edit29 handles.edit30 handles.edit31 handles.edit32 handles.edit33 handles.edit34 handles.edit35 handles.edit36 handles.edit37 handles.edit38 handles.edit39 handles.edit40 handles.edit41 handles.edit42 handles.edit43 handles.edit44 handles.edit45 handles.edit46 handles.edit47 handles.edit48 handles.edit49 handles.edit50 handles.edit51 handles.edit52 handles.edit53 handles.edit54 handles.edit55 handles.edit56 handles.edit57 handles.edit58 handles.edit59 handles.edit60 handles.edit61 handles.edit62 handles.edit63 handles.edit64 handles.edit65 handles.edit66 handles.edit67 handles.edit68 handles.edit69 handles.edit70 handles.edit71 handles.edit72 handles.edit73 handles.edit74 handles.edit75 handles.edit76 handles.edit77 handles.edit78 handles.edit79 handles.edit80 handles.edit81 handles.edit82 handles.edit83 handles.edit84];
    
    
    % SET TIME WINDOW in ms, except 0 which mean no margin
    set(handles.slider1,'Value',200)
    set(handles.slider1,'Min',0)
    set(handles.slider1,'Max',300)
    
    % DISPLAY THE AXES
    for i=1:84
        if (handles.page-1)*84+i < length(handles.syltempType)
            set(handles.edit(i),'String',handles.syltempType((handles.page-1)*84+i))
            set(handles.edit(i), 'ButtonDownFcn',@cb_remove_file)
            set(handles.edit(i),'BackgroundColor',[1 1 1])
        else
            set(handles.edit(i),'String',0)
            set(handles.edit(i),'BackgroundColor',[1 1 1])
        end
        handles = updateAudAxes(handles,(handles.page-1)*84+i,handles.axes(i));
    end
    
end

guidata(hObject, handles);

% --- Executes on button press in push_load.
function push_load_Callback(hObject, eventdata, handles)

directory_name = uigetdir();
if ~ischar(directory_name)
    return
else
    handles.dirname = directory_name;
    cd(handles.dirname);
    handles.filelist = dir('*.wav');
end

guidata(hObject, handles);

function clickonbox(hObject)
handles = guidata(hObject);
set(hObject,'BackgroundColor',[1 0.6 0.6])
guidata(hObject, handles);

function handles = updateAudAxes(handles,IndxFile,axes_spec)
axes(axes_spec);
if IndxFile > length(handles.syltempType)
    cla;
    xlabel('');
    ylabel('');
else
    if get(handles.slider1,'Value')~=0
    winsize=get(handles.slider1,'Value')*44.15;
    margin=floor(winsize-(handles.syltempEnd(IndxFile)-handles.syltempStart(IndxFile))/2);
    stri=floor(handles.syltempStart(IndxFile)-margin);
    endi=floor(handles.syltempEnd(IndxFile)+margin);
    else
    margin=800;
    stri=floor(handles.syltempStart(IndxFile)-margin);
    endi=floor(handles.syltempEnd(IndxFile)+margin);
    end
        
    try
        audio=handles.audio{handles.syltempIndxKey(IndxFile)}(stri:endi);
    catch err
        if stri<0            
        audio=handles.audio{handles.syltempIndxKey(IndxFile)};
        audio=vertcat(zeros(-stri+1,1), audio(1:endi));
        else
        audio=handles.audio{handles.syltempIndxKey(IndxFile)};
        audio=vertcat(audio(stri:end), zeros((endi-length(audio)),1));
        end
    end
    displaySpecgramQuickalex(audio, 44150,[0,10000],[],0);
    if get(handles.slider1,'Value')~=0
    minimargin=300;
    x = [(margin-minimargin)/44150, (margin-minimargin)/44150, (margin+handles.syltempEnd(IndxFile)-handles.syltempStart(IndxFile)+minimargin)/44150, (margin+handles.syltempEnd(IndxFile)-handles.syltempStart(IndxFile)+minimargin)/44150, (margin-minimargin)/44150];
    lims = ylim;
    y = [lims(1)-2, lims(2)+2, lims(2)+2, lims(1)-2, lims(1)-2];
    patch(x,y,[1 0 0], 'EdgeColor', [1 0 0], 'FaceAlpha', 0, 'LineWidth', 1.1);
    end
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca, 'ButtonDownFcn',@cb_display_date)
    xlabel('');
    ylabel('');
end



% --- Executes on selection change in popupmenu.
function popupmenu_Callback(hObject, eventdata, handles)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu
if ~isempty(handles.syltempType)
    for i=1:84
        if str2double(get(handles.edit(i),'String')) ~= 0
        handles.syltempType((handles.page-1)*84+i)= str2double(get(handles.edit(i),'String'));
        end
    end
    handles.sylType(handles.syltempIndx)=handles.syltempType(1:length(handles.sylType(handles.syltempIndx)));
end
    handles.page=1;
    contents = cellstr(get(handles.popupmenu,'String'));
    handles.syltempvalue=str2double(contents{get(handles.popupmenu,'Value')});
    handles.syltempIndx=find(handles.sylType==handles.syltempvalue);
    handles.syltempBin=handles.sylType==handles.syltempvalue;
    handles.syltempStart=handles.sylStart(handles.syltempBin);
    handles.syltempEnd=handles.sylEnd(handles.syltempBin);
    handles.syltempType=handles.sylType(handles.syltempBin);
    handles.syltempIndxKey=handles.IndxKey(handles.syltempBin);
    handles.temp_max_page=floor(length(handles.syltempType)/84);
    set(handles.popupmenu_page,'String',num2str((1:(handles.temp_max_page+1))'));
    set(handles.popupmenu_page,'Value',1);

contents = cellstr(get(handles.popupmenu,'String'));
sylnum=str2double(contents{get(handles.popupmenu,'Value')});

for i=1:84
    handles = updateAudAxes(handles,(handles.page-1)*84+i,handles.axes(i));
    if (handles.page-1)*84+i <= length(handles.syltempType)
        set(handles.edit(i),'String',handles.syltempType((handles.page-1)*84+i))
        if handles.syltempType((handles.page-1)*84+i) ~= sylnum
            set(handles.edit(i),'BackgroundColor',[1 0.6 0.6])
        else
                set(handles.edit(i),'BackgroundColor',[1 1 1])
        end     
    else
        set(handles.edit(i),'String',0)
        set(handles.edit(i),'BackgroundColor',[1 1 1])
    end
end

guidata(hObject, handles);



function popupmenu_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_next.
function push_next_Callback(hObject, eventdata, handles)

for i=1:84
    handles.syltempType((handles.page-1)*84+i)= str2double(get(handles.edit(i),'String'));
end
handles.sylType(handles.syltempIndx)=handles.syltempType(1:length(handles.sylType(handles.syltempIndx)));


handles.page=handles.page+1;
contents = cellstr(get(handles.popupmenu_page,'String'));
set(handles.popupmenu_page,'Value',str2double(contents{get(handles.popupmenu_page,'Value')})+1);
handles.datedisplay=cell(1,84);


contents = cellstr(get(handles.popupmenu,'String'));
sylnum=str2double(contents{get(handles.popupmenu,'Value')});

for i=1:84
    handles = updateAudAxes(handles,(handles.page-1)*84+i,handles.axes(i));
    if (handles.page-1)*84+i <= length(handles.syltempType)
        set(handles.edit(i),'String',handles.syltempType((handles.page-1)*84+i))
        if handles.syltempType((handles.page-1)*84+i) ~= sylnum
            set(handles.edit(i),'BackgroundColor',[1 0.6 0.6])
        else
                set(handles.edit(i),'BackgroundColor',[1 1 1])
        end        
    else
        set(handles.edit(i),'String',0)
        set(handles.edit(i),'BackgroundColor',[1 1 1])
    end
end


guidata(hObject, handles);


% --- Executes on button press in push_previous.
function push_previous_Callback(hObject, eventdata, handles)
if handles.page ~=1
    
    for i=1:84
        if str2double(get(handles.edit(i),'String')) ~= 0
        handles.syltempType((handles.page-1)*84+i)= str2double(get(handles.edit(i),'String'));
        end
    end
    handles.sylType(handles.syltempIndx)=handles.syltempType(1:length(handles.sylType(handles.syltempIndx)));
    
    
    handles.page=handles.page-1;
    contents = cellstr(get(handles.popupmenu_page,'String'));
    set(handles.popupmenu_page,'Value',str2double(contents{get(handles.popupmenu_page,'Value')})-1);
     
    
contents = cellstr(get(handles.popupmenu,'String'));
sylnum=str2double(contents{get(handles.popupmenu,'Value')});

    
    for i=1:84
        handles = updateAudAxes(handles,(handles.page-1)*84+i,handles.axes(i));
        if (handles.page-1)*84+i <= length(handles.syltempType)
            set(handles.edit(i),'String',handles.syltempType((handles.page-1)*84+i))
            if handles.syltempType((handles.page-1)*84+i) ~= sylnum
                set(handles.edit(i),'BackgroundColor',[1 0.6 0.6])
            else
                set(handles.edit(i),'BackgroundColor',[1 1 1])
            end      
        else
            set(handles.edit(i),'String',0)
        end
    end
end
guidata(hObject, handles);



% --- Executes on button press in push_viewsylpattern.
function push_viewsylpattern_Callback(hObject, eventdata, handles)
if isfield(handles,'templates')
    ansButton = questdlg('A template file is already loaded. Keep this one or load another?','Template Already Loaded','Keep & View','Load Another','Keep & View');
    if strcmp(ansButton,'Load Another')
        uiload;
        handles.templates = templates;
    end
else
    uiload;
    handles.templates = templates;
end
handles = updateTemplates(handles);
guidata(hObject, handles);



function handles = updateTemplates(handles)
h = figure(6969);
position = get(h,'Position');
handles.axesTemplates = axes;
cla(handles.axesTemplates);

if(isfield(handles, 'templates'))
    stdfs = 44150;
    allwavs = [];
    edges = [0];
    for(nWav = 1:length(handles.templates.wavs))
        fs = handles.templates.wavs(nWav).fs;
        wav = handles.templates.wavs(nWav).wav;
        label = handles.templates.wavs(nWav).segType;
        wav = resample(wav, stdfs, fs);
        allwavs = [allwavs; wav'];
        edges = [edges, edges(end) + (length(wav) - 1)/stdfs];
    end
    
    %zero pad to 1 sec in length.
    if(edges(end) < 1)
        allwavs(end+1:stdfs) = 0;
    end
    
    %display the templates
    axes(handles.axesTemplates);
    displaySpecgramQuick(allwavs, stdfs, [0,10000]);
    
    %seperate with lines and label.
    for(nWav = 1:length(edges)-1)
        l = line([edges(nWav),edges(nWav)], ylim);
        set(l,'Color','k');
        set(l,'LineWidth', 2);
        mid = mean(edges([nWav,nWav+1]));
        if(handles.templates.wavs(nWav).segType ~= -1)
            t = text(mid, mean(ylim), num2str(handles.templates.wavs(nWav).segType), 'FontWeight', 'bold', 'FontSize', 20);
        end
    end
    
    set(handles.axesTemplates, 'ButtonDownFcn', @cb_templates_click);
    position(1,3) = 1024; position(1,4) = 180;
    set(h,'Position',position);
end


% --- Executes on button press in push_delRow1.
function push_delRow1_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
set(handles.edit(1:12),'String',get(handles.edit_rightClick,'String'))
set(handles.edit(1:12),'BackgroundColor',[1 0.6 0.6])
guidata(hObject, handles);

% --- Executes on button press in push_delRow2.
function push_delRow2_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
set(handles.edit(13:24),'String',get(handles.edit_rightClick,'String'))
set(handles.edit(13:24),'BackgroundColor',[1 0.6 0.6])
guidata(hObject, handles);

% --- Executes on button press in push_delRow3.
function push_delRow3_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
set(handles.edit(25:36),'String',get(handles.edit_rightClick,'String'))
set(handles.edit(25:36),'BackgroundColor',[1 0.6 0.6])
guidata(hObject, handles);

% --- Executes on button press in push_delRow4.
function push_delRow4_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
set(handles.edit(37:48),'String',get(handles.edit_rightClick,'String'))
set(handles.edit(37:48),'BackgroundColor',[1 0.6 0.6])
guidata(hObject, handles);

% --- Executes on button press in push_delRow5.
function push_delRow5_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
set(handles.edit(49:60),'String',get(handles.edit_rightClick,'String'))
set(handles.edit(49:60),'BackgroundColor',[1 0.6 0.6])
guidata(hObject, handles);

% --- Executes on button press in push_delRow6.
function push_delRow6_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
set(handles.edit(61:72),'String',get(handles.edit_rightClick,'String'))
set(handles.edit(61:72),'BackgroundColor',[1 0.6 0.6])
guidata(hObject, handles);

% --- Executes on button press in push_delRow7.
function push_delRow7_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
set(handles.edit(72:84),'String',get(handles.edit_rightClick,'String'))
set(handles.edit(72:84),'BackgroundColor',[1 0.6 0.6])
guidata(hObject, handles);




% --- Executes on selection change in popupmenu_page.
function popupmenu_page_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_page contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_page

for i=1:84
    if str2double(get(handles.edit(i),'String')) ~= 0
        handles.syltempType((handles.page-1)*84+i)= str2double(get(handles.edit(i),'String'));
    end
end
handles.sylType(handles.syltempIndx)=handles.syltempType(1:length(handles.sylType(handles.syltempIndx)));

contents = cellstr(get(handles.popupmenu_page,'String'));
handles.page=str2double(contents{get(handles.popupmenu_page,'Value')});

contentssyl = cellstr(get(handles.popupmenu,'String'));
sylnum=str2double(contentssyl{get(handles.popupmenu,'Value')});

for i=1:84
    handles = updateAudAxes(handles,(handles.page-1)*84+i,handles.axes(i));
    if (handles.page-1)*84+i <= length(handles.syltempType)
        set(handles.edit(i),'String',handles.syltempType((handles.page-1)*84+i))
        if handles.syltempType((handles.page-1)*84+i) ~= sylnum
            set(handles.edit(i),'BackgroundColor',[1 0.6 0.6])
        else
                set(handles.edit(i),'BackgroundColor',[1 1 1])
        end  
    else
        set(handles.edit(i),'String',0)
    end
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_page_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_rightClick_Callback(hObject, eventdata, handles)
function edit_rightClick_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_sylclus.
function push_sylclus_Callback(hObject, eventdata, handles)
% hObject    handle to push_sylclus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% SAVE DATA
if ~isempty(handles.syltempType)
    for i=1:84
        if str2double(get(handles.edit(i),'String')) ~= 0
        handles.syltempType((handles.page-1)*84+i)= str2double(get(handles.edit(i),'String'));
        end
    end
    handles.sylType(handles.syltempIndx)=handles.syltempType(1:length(handles.sylType(handles.syltempIndx)));
end


% LAUNCH SYLCLUS
input_sylclus.audio=handles.audio;
input_sylclus.sylStart=handles.sylStart;
input_sylclus.sylEnd=handles.sylEnd;
input_sylclus.sylType=handles.sylType;
input_sylclus.IndxKey=handles.IndxKey;
contents = cellstr(get(handles.popupmenu,'String'));
input_sylclus.sylvalue=str2double(contents{get(handles.popupmenu,'Value')});

output_sylclus = SylClus(input_sylclus);

handles.sylStart=output_sylclus.sylStart ;
handles.sylEnd=output_sylclus.sylEnd;
handles.sylType=output_sylclus.sylType;
handles.IndxKey=output_sylclus.IndxKey;
handles.page=1;

%DISPLAY NEW RESULT
contentssyl = cellstr(get(handles.popupmenu,'String'));
sylnum=str2double(contentssyl{get(handles.popupmenu,'Value')});
    contents = cellstr(get(handles.popupmenu,'String'));
    handles.syltempvalue=str2double(contents{get(handles.popupmenu,'Value')});
    handles.syltempIndx=find(handles.sylType==handles.syltempvalue);
    handles.syltempBin=handles.sylType==handles.syltempvalue;
    handles.syltempStart=handles.sylStart(handles.syltempBin);
    handles.syltempEnd=handles.sylEnd(handles.syltempBin);
    handles.syltempType=handles.sylType(handles.syltempBin);
    handles.syltempIndxKey=handles.IndxKey(handles.syltempBin);
    handles.temp_max_page=floor(length(handles.syltempType)/84);
    set(handles.popupmenu_page,'String',num2str((1:(handles.temp_max_page+1))'));
    set(handles.popupmenu_page,'Value',1);
    for i=1:84
    handles = updateAudAxes(handles,(handles.page-1)*84+i,handles.axes(i));
    if (handles.page-1)*84+i <= length(handles.syltempType)
        set(handles.edit(i),'String',handles.syltempType((handles.page-1)*84+i))
        if handles.syltempType((handles.page-1)*84+i) ~= sylnum
            set(handles.edit(i),'BackgroundColor',[1 0.6 0.6])
        else
                set(handles.edit(i),'BackgroundColor',[1 1 1])
        end     
    else
        set(handles.edit(i),'String',0)
        set(handles.edit(i),'BackgroundColor',[1 1 1])
    end
end


guidata(hObject, handles);

































function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clickonbox(hObject); guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit28 as text
%        str2double(get(hObject,'String')) returns contents of edit28 as a double


% --- Executes during object creation, after setting all properties.
function edit28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit29_Callback(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit29 as text
%        str2double(get(hObject,'String')) returns contents of edit29 as a double


% --- Executes during object creation, after setting all properties.
function edit29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit30_Callback(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit30 as text
%        str2double(get(hObject,'String')) returns contents of edit30 as a double


% --- Executes during object creation, after setting all properties.
function edit30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit31_Callback(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit31 as text
%        str2double(get(hObject,'String')) returns contents of edit31 as a double


% --- Executes during object creation, after setting all properties.
function edit31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit32_Callback(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit32 as text
%        str2double(get(hObject,'String')) returns contents of edit32 as a double


% --- Executes during object creation, after setting all properties.
function edit32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit33_Callback(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit33 as text
%        str2double(get(hObject,'String')) returns contents of edit33 as a double


% --- Executes during object creation, after setting all properties.
function edit33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit34_Callback(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit34 as text
%        str2double(get(hObject,'String')) returns contents of edit34 as a double


% --- Executes during object creation, after setting all properties.
function edit34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit35_Callback(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit35 as text
%        str2double(get(hObject,'String')) returns contents of edit35 as a double


% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit36_Callback(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit36 as text
%        str2double(get(hObject,'String')) returns contents of edit36 as a double


% --- Executes during object creation, after setting all properties.
function edit36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit37_Callback(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit37 as text
%        str2double(get(hObject,'String')) returns contents of edit37 as a double


% --- Executes during object creation, after setting all properties.
function edit37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit38_Callback(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit38 as text
%        str2double(get(hObject,'String')) returns contents of edit38 as a double


% --- Executes during object creation, after setting all properties.
function edit38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit39_Callback(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit39 as text
%        str2double(get(hObject,'String')) returns contents of edit39 as a double


% --- Executes during object creation, after setting all properties.
function edit39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit40_Callback(hObject, eventdata, handles)
% hObject    handle to edit40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit40 as text
%        str2double(get(hObject,'String')) returns contents of edit40 as a double


% --- Executes during object creation, after setting all properties.
function edit40_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit41_Callback(hObject, eventdata, handles)
% hObject    handle to edit41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit41 as text
%        str2double(get(hObject,'String')) returns contents of edit41 as a double


% --- Executes during object creation, after setting all properties.
function edit41_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit42_Callback(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit42 as text
%        str2double(get(hObject,'String')) returns contents of edit42 as a double


% --- Executes during object creation, after setting all properties.
function edit42_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit43_Callback(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit43 as text
%        str2double(get(hObject,'String')) returns contents of edit43 as a double


% --- Executes during object creation, after setting all properties.
function edit43_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit44_Callback(hObject, eventdata, handles)
% hObject    handle to edit44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit44 as text
%        str2double(get(hObject,'String')) returns contents of edit44 as a double


% --- Executes during object creation, after setting all properties.
function edit44_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit45_Callback(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit45 as text
%        str2double(get(hObject,'String')) returns contents of edit45 as a double


% --- Executes during object creation, after setting all properties.
function edit45_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit46_Callback(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit46 as text
%        str2double(get(hObject,'String')) returns contents of edit46 as a double


% --- Executes during object creation, after setting all properties.
function edit46_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit47_Callback(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit47 as text
%        str2double(get(hObject,'String')) returns contents of edit47 as a double


% --- Executes during object creation, after setting all properties.
function edit47_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit48_Callback(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit48 as text
%        str2double(get(hObject,'String')) returns contents of edit48 as a double


% --- Executes during object creation, after setting all properties.
function edit48_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit49_Callback(hObject, eventdata, handles)
% hObject    handle to edit49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit49 as text
%        str2double(get(hObject,'String')) returns contents of edit49 as a double


% --- Executes during object creation, after setting all properties.
function edit49_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit50_Callback(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit50 as text
%        str2double(get(hObject,'String')) returns contents of edit50 as a double


% --- Executes during object creation, after setting all properties.
function edit50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit51_Callback(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit51 as text
%        str2double(get(hObject,'String')) returns contents of edit51 as a double


% --- Executes during object creation, after setting all properties.
function edit51_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit52_Callback(hObject, eventdata, handles)
% hObject    handle to edit52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit52 as text
%        str2double(get(hObject,'String')) returns contents of edit52 as a double


% --- Executes during object creation, after setting all properties.
function edit52_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit53_Callback(hObject, eventdata, handles)
% hObject    handle to edit53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit53 as text
%        str2double(get(hObject,'String')) returns contents of edit53 as a double


% --- Executes during object creation, after setting all properties.
function edit53_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit54_Callback(hObject, eventdata, handles)
% hObject    handle to edit54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit54 as text
%        str2double(get(hObject,'String')) returns contents of edit54 as a double


% --- Executes during object creation, after setting all properties.
function edit54_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit55_Callback(hObject, eventdata, handles)
% hObject    handle to edit55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit55 as text
%        str2double(get(hObject,'String')) returns contents of edit55 as a double


% --- Executes during object creation, after setting all properties.
function edit55_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit56_Callback(hObject, eventdata, handles)
% hObject    handle to edit56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit56 as text
%        str2double(get(hObject,'String')) returns contents of edit56 as a double


% --- Executes during object creation, after setting all properties.
function edit56_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit57_Callback(hObject, eventdata, handles)
% hObject    handle to edit57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit57 as text
%        str2double(get(hObject,'String')) returns contents of edit57 as a double


% --- Executes during object creation, after setting all properties.
function edit57_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit58_Callback(hObject, eventdata, handles)
% hObject    handle to edit58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit58 as text
%        str2double(get(hObject,'String')) returns contents of edit58 as a double


% --- Executes during object creation, after setting all properties.
function edit58_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit59_Callback(hObject, eventdata, handles)
% hObject    handle to edit59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit59 as text
%        str2double(get(hObject,'String')) returns contents of edit59 as a double


% --- Executes during object creation, after setting all properties.
function edit59_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit60_Callback(hObject, eventdata, handles)
% hObject    handle to edit60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit60 as text
%        str2double(get(hObject,'String')) returns contents of edit60 as a double


% --- Executes during object creation, after setting all properties.
function edit60_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit61_Callback(hObject, eventdata, handles)
% hObject    handle to edit61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit61 as text
%        str2double(get(hObject,'String')) returns contents of edit61 as a double


% --- Executes during object creation, after setting all properties.
function edit61_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit62_Callback(hObject, eventdata, handles)
% hObject    handle to edit62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit62 as text
%        str2double(get(hObject,'String')) returns contents of edit62 as a double


% --- Executes during object creation, after setting all properties.
function edit62_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit63_Callback(hObject, eventdata, handles)
% hObject    handle to edit63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit63 as text
%        str2double(get(hObject,'String')) returns contents of edit63 as a double


% --- Executes during object creation, after setting all properties.
function edit63_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit64_Callback(hObject, eventdata, handles)
% hObject    handle to edit64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit64 as text
%        str2double(get(hObject,'String')) returns contents of edit64 as a double


% --- Executes during object creation, after setting all properties.
function edit64_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit65_Callback(hObject, eventdata, handles)
% hObject    handle to edit65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit65 as text
%        str2double(get(hObject,'String')) returns contents of edit65 as a double


% --- Executes during object creation, after setting all properties.
function edit65_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit66_Callback(hObject, eventdata, handles)
% hObject    handle to edit66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit66 as text
%        str2double(get(hObject,'String')) returns contents of edit66 as a double


% --- Executes during object creation, after setting all properties.
function edit66_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit67_Callback(hObject, eventdata, handles)
% hObject    handle to edit67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit67 as text
%        str2double(get(hObject,'String')) returns contents of edit67 as a double


% --- Executes during object creation, after setting all properties.
function edit67_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit68_Callback(hObject, eventdata, handles)
% hObject    handle to edit68 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit68 as text
%        str2double(get(hObject,'String')) returns contents of edit68 as a double


% --- Executes during object creation, after setting all properties.
function edit68_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit68 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit69_Callback(hObject, eventdata, handles)
% hObject    handle to edit69 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit69 as text
%        str2double(get(hObject,'String')) returns contents of edit69 as a double


% --- Executes during object creation, after setting all properties.
function edit69_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit69 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit70_Callback(hObject, eventdata, handles)
% hObject    handle to edit70 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit70 as text
%        str2double(get(hObject,'String')) returns contents of edit70 as a double


% --- Executes during object creation, after setting all properties.
function edit70_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit70 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit71_Callback(hObject, eventdata, handles)
% hObject    handle to edit71 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit71 as text
%        str2double(get(hObject,'String')) returns contents of edit71 as a double


% --- Executes during object creation, after setting all properties.
function edit71_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit71 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit72_Callback(hObject, eventdata, handles)
% hObject    handle to edit72 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit72 as text
%        str2double(get(hObject,'String')) returns contents of edit72 as a double


% --- Executes during object creation, after setting all properties.
function edit72_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit72 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit73_Callback(hObject, eventdata, handles)
% hObject    handle to edit73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit73 as text
%        str2double(get(hObject,'String')) returns contents of edit73 as a double


% --- Executes during object creation, after setting all properties.
function edit73_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit74_Callback(hObject, eventdata, handles)
% hObject    handle to edit74 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit74 as text
%        str2double(get(hObject,'String')) returns contents of edit74 as a double


% --- Executes during object creation, after setting all properties.
function edit74_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit74 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit75_Callback(hObject, eventdata, handles)
% hObject    handle to edit75 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit75 as text
%        str2double(get(hObject,'String')) returns contents of edit75 as a double


% --- Executes during object creation, after setting all properties.
function edit75_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit75 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit76_Callback(hObject, eventdata, handles)
% hObject    handle to edit76 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit76 as text
%        str2double(get(hObject,'String')) returns contents of edit76 as a double


% --- Executes during object creation, after setting all properties.
function edit76_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit76 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit77_Callback(hObject, eventdata, handles)
% hObject    handle to edit77 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit77 as text
%        str2double(get(hObject,'String')) returns contents of edit77 as a double


% --- Executes during object creation, after setting all properties.
function edit77_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit77 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit78_Callback(hObject, eventdata, handles)
% hObject    handle to edit78 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit78 as text
%        str2double(get(hObject,'String')) returns contents of edit78 as a double


% --- Executes during object creation, after setting all properties.
function edit78_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit78 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit79_Callback(hObject, eventdata, handles)
% hObject    handle to edit79 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit79 as text
%        str2double(get(hObject,'String')) returns contents of edit79 as a double


% --- Executes during object creation, after setting all properties.
function edit79_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit79 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit80_Callback(hObject, eventdata, handles)
% hObject    handle to edit80 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit80 as text
%        str2double(get(hObject,'String')) returns contents of edit80 as a double


% --- Executes during object creation, after setting all properties.
function edit80_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit80 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit81_Callback(hObject, eventdata, handles)
% hObject    handle to edit81 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit81 as text
%        str2double(get(hObject,'String')) returns contents of edit81 as a double


% --- Executes during object creation, after setting all properties.
function edit81_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit81 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit82_Callback(hObject, eventdata, handles)
% hObject    handle to edit82 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit82 as text
%        str2double(get(hObject,'String')) returns contents of edit82 as a double


% --- Executes during object creation, after setting all properties.
function edit82_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit82 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit83_Callback(hObject, eventdata, handles)
% hObject    handle to edit83 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit83 as text
%        str2double(get(hObject,'String')) returns contents of edit83 as a double


% --- Executes during object creation, after setting all properties.
function edit83_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit83 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit84_Callback(hObject, eventdata, handles)
% hObject    handle to edit84 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
clickonbox(hObject); guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit84 as text
%        str2double(get(hObject,'String')) returns contents of edit84 as a double


% --- Executes during object creation, after setting all properties.
function edit84_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit84 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
