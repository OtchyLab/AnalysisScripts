function varargout = SylClus(varargin)
% SYLCLUS MATLAB code for SylClus.fig
%      SYLCLUS, by itself, creates a new SYLCLUS or raises the existing
%      singleton*.
%
%      H = SYLCLUS returns the handle to a new SYLCLUS or the handle to
%      the existing singleton*.
%
%      SYLCLUS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SYLCLUS.M with the given input arguments.
%
%      SYLCLUS('Property','Value',...) creates a new SYLCLUS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SylClus_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SylClus_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SylClus

% Last Modified by GUIDE v2.5 18-Mar-2014 12:54:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SylClus_OpeningFcn, ...
                   'gui_OutputFcn',  @SylClus_OutputFcn, ...
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


% --- Executes just before SylClus is made visible.
function SylClus_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SylClus (see VARARGIN)

% PASTE IN HANDLES
handles.audio = varargin{1}.audio;
handles.sylStart = varargin{1}.sylStart;
handles.sylEnd = varargin{1}.sylEnd;
handles.sylType = varargin{1}.sylType;
handles.IndxKey = varargin{1}.IndxKey;
handles.sylvalue = varargin{1}.sylvalue;

% PARAMETER FOR THE DISPLAY

handles.syltempIndx=find(handles.sylType==handles.sylvalue);
handles.syltempBin=handles.sylType==handles.sylvalue;
handles.syltempStart=handles.sylStart(handles.syltempBin);
handles.syltempEnd=handles.sylEnd(handles.syltempBin);
handles.syltempType=handles.sylType(handles.syltempBin);
handles.syltempIndxKey=handles.IndxKey(handles.syltempBin);


%CLEAR PLOT
axes(handles.axes_syl1)
cla
axes(handles.axes_syl2)
cla
axes(handles.axes_syl3)
cla
axes(handles.axes_syl4)
cla

% SPECTRO INTERPOL
SpecSyl=cell(length(handles.syltempStart),1);
margin= floor(0.01*44150);
time=zeros(length(SpecSyl),1);

startmargin=floor(handles.syltempStart-margin);
startmargin(startmargin<=0)=1;
endmargin=floor(handles.syltempEnd+margin);

%CALCUL SPECTRO FOR EACH FILE
parfor i=1:length(handles.syltempStart)
    try
        SpecSyl{i}=handles.audio{handles.syltempIndxKey(i)}(startmargin(i):endmargin(i));
    catch 
        SpecSyl{i}=handles.audio{handles.syltempIndxKey(i)}(startmargin(i):end);
    end
    [~,~,~,SpecSyl{i}] = spectrogram((SpecSyl{i}/(sqrt(mean(SpecSyl{i}.^2)))),220,220-44,512,44150);
    strt = 4; stp = 94; %These bins correspond to 300-8000Hz
    SpecSyl{i} = abs(SpecSyl{i}(strt:stp,:));
    time(i)=size(SpecSyl{i},2);
end

% INTERPOL SPECTRO
freq=size(SpecSyl{1},1);
tim=floor(mean(time));

SpecSylinter=zeros(size(SpecSyl{1},1),tim ,length(SpecSyl));


parfor i=1:length(SpecSyl)
    for j=1:freq
    SpecSylinter(j,:,i)=interp1(SpecSyl{i}(j,:),linspace(1,size(SpecSyl{i},2),tim));
    end
end

handles.SpecSylinter=SpecSylinter;

% CRITERION TO DISPLAY
crit_bin=(handles.syltempEnd-handles.syltempStart)>0;
MeanSpec=log10(mean(SpecSylinter(:,:,crit_bin),3));

axes(handles.axes_syl1)
imagesc(MeanSpec)
axis xy


handles.limits=sort([get(handles.slider_bound1,'Value') get(handles.slider_bound2,'Value') get(handles.slider_bound3,'Value')]);
duration=(handles.syltempEnd-handles.syltempStart)/44150;

% SET SLIDER AND EDIT BOX
%%%%%%%%%%%%%%%%%%%%
set(handles.slider_classnumber,'Value',100)
set(handles.slider_classnumber,'Min',1)
set(handles.slider_classnumber,'Max',800)

set(handles.slider_bound1,'Value',0)
set(handles.slider_bound1,'Min',0)
set(handles.slider_bound1,'Max',max(duration)+0.01*max(duration))
set(handles.slider_bound1, 'SliderStep', [0.0005 0.001])

set(handles.slider_bound2,'Value',0)
set(handles.slider_bound2,'Min',0)
set(handles.slider_bound2,'Max',max(duration)+0.01*max(duration))
set(handles.slider_bound2, 'SliderStep', [0.0005 0.001])

set(handles.slider_bound3,'Value',0)
set(handles.slider_bound3,'Min',0)
set(handles.slider_bound3,'Max',max(duration)+0.01*max(duration))
set(handles.slider_bound3, 'SliderStep', [0.0005 0.001])

set(handles.edit_bound1,'String',get(handles.slider_bound1,'Value'))
set(handles.edit_bound1,'BackgroundColor', [1 0.6 0.6])
set(handles.edit_bound2,'String',get(handles.slider_bound2,'Value'))
set(handles.edit_bound2,'BackgroundColor', [0.6 1 0.6])
set(handles.edit_bound3,'String',get(handles.slider_bound3,'Value'))
set(handles.edit_bound3,'BackgroundColor', [0.6 0.6 0.6])

set(handles.edit_rename1,'String',handles.sylvalue)
set(handles.edit_rename2,'String',handles.sylvalue)
set(handles.edit_rename3,'String',handles.sylvalue)
set(handles.edit_rename4,'String',handles.sylvalue)

set(handles.edit_split_left1,'String',handles.sylvalue)
set(handles.edit_split_left2,'String',handles.sylvalue)
set(handles.edit_split_left3,'String',handles.sylvalue)
set(handles.edit_split_left4,'String',handles.sylvalue)

set(handles.edit_split_right1,'String',handles.sylvalue)
set(handles.edit_split_right2,'String',handles.sylvalue)
set(handles.edit_split_right3,'String',handles.sylvalue)
set(handles.edit_split_right4,'String',handles.sylvalue)

set(handles.slider_split1,'Value',0)
set(handles.slider_split1,'Min',0)
set(handles.slider_split1,'Max',tim)
set(handles.slider_split1, 'SliderStep', [0.1 1]) 

set(handles.slider_split2,'Value',0)
set(handles.slider_split2,'Min',0)
set(handles.slider_split2,'Max',tim)
set(handles.slider_split2, 'SliderStep', [0.1 1])

set(handles.slider_split3,'Value',0)
set(handles.slider_split3,'Min',0)
set(handles.slider_split3,'Max',tim)
set(handles.slider_split3, 'SliderStep', [0.1 1]) 

set(handles.slider_split4,'Value',0)
set(handles.slider_split4,'Min',0)
set(handles.slider_split4,'Max',tim)
set(handles.slider_split4, 'SliderStep', [0.1 1]) 

% DISPLAY HISTO
handles=refresh_hist(handles);



% Choose default command line output for SylClus
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SylClus wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SylClus_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

duration=(handles.syltempEnd-handles.syltempStart)/44150;
handles.sylType(handles.syltempIndx(duration<=handles.limits(1)))=str2double(get(handles.edit_rename1,'String'));
handles.sylType(handles.syltempIndx((duration<=handles.limits(2)) & (duration>handles.limits(1))))=str2double(get(handles.edit_rename2,'String'));
handles.sylType(handles.syltempIndx((duration<=handles.limits(3)) & (duration>handles.limits(2))))=str2double(get(handles.edit_rename3,'String'));
handles.sylType(handles.syltempIndx((duration>handles.limits(3))))=str2double(get(handles.edit_rename4,'String'));

output_sylclus.sylStart=handles.sylStart ;
output_sylclus.sylEnd=handles.sylEnd ;
output_sylclus.sylType=handles.sylType ;
output_sylclus.IndxKey=handles.IndxKey ;

varargout{1}=output_sylclus;
delete(handles.figure1);


% --- Executes on slider movement.
function slider_classnumber_Callback(hObject, eventdata, handles)
% hObject    handle to slider_classnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles=refresh_hist(handles);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_classnumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_classnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_bound1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bound1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bound1 as text
%        str2double(get(hObject,'String')) returns contents of edit_bound1 as a double

handles=refresh_hist(handles);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_bound1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bound1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_bound1_Callback(hObject, eventdata, handles)
% hObject    handle to slider_bound1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(handles.edit_bound1,'String',get(handles.slider_bound1,'Value'))
handles=refresh_hist(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_bound1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_bound1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_bound2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bound2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bound2 as text
%        str2double(get(hObject,'String')) returns contents of edit_bound2 as a double
handles=refresh_hist(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_bound2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bound2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_bound2_Callback(hObject, eventdata, handles)
% hObject    handle to slider_bound2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(handles.edit_bound2,'String',get(handles.slider_bound2,'Value'))
handles=refresh_hist(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_bound2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_bound2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_bound3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bound3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bound3 as text
%        str2double(get(hObject,'String')) returns contents of edit_bound3 as a double
handles=refresh_hist(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_bound3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bound3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_bound3_Callback(hObject, eventdata, handles)
% hObject    handle to slider_bound3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.edit_bound3,'String',get(handles.slider_bound3,'Value'))
handles=refresh_hist(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_bound3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_bound3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in push_refresh.
function push_refresh_Callback(hObject, eventdata, handles)
% hObject    handle to push_refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.syltempIndx=find(handles.sylType==handles.sylvalue);
handles.syltempBin=handles.sylType==handles.sylvalue;
handles.syltempStart=handles.sylStart(handles.syltempBin);
handles.syltempEnd=handles.sylEnd(handles.syltempBin);
handles.syltempType=handles.sylType(handles.syltempBin);
handles.syltempIndxKey=handles.IndxKey(handles.syltempBin);


handles.limits=sort([get(handles.slider_bound1,'Value') get(handles.slider_bound2,'Value') get(handles.slider_bound3,'Value')]);
duration=(handles.syltempEnd-handles.syltempStart)/44150;


handles.MeanSpec=[];
crit_bin=duration>handles.limits(3);
handles.MeanSpec(:,:,4)=log10(mean(handles.SpecSylinter(:,:,crit_bin),3));
axes(handles.axes_syl4)
imagesc(handles.MeanSpec(:,:,4))
axis xy

crit_bin=(duration<=handles.limits(3)) & (duration>handles.limits(2));
handles.MeanSpec(:,:,3)=log10(mean(handles.SpecSylinter(:,:,crit_bin),3));
axes(handles.axes_syl3)
imagesc(handles.MeanSpec(:,:,3))
axis xy

crit_bin=(duration<=handles.limits(2)) & (duration>handles.limits(1));
handles.MeanSpec(:,:,2)=log10(mean(handles.SpecSylinter(:,:,crit_bin),3));
axes(handles.axes_syl2)
imagesc(handles.MeanSpec(:,:,2))
axis xy

crit_bin=duration<=handles.limits(1);
handles.MeanSpec(:,:,1)=log10(mean(handles.SpecSylinter(:,:,crit_bin),3));
axes(handles.axes_syl1)
imagesc(handles.MeanSpec(:,:,1))
axis xy

set(handles.edit_rename1,'String',handles.sylvalue)
set(handles.edit_rename2,'String',handles.sylvalue)
set(handles.edit_rename3,'String',handles.sylvalue)
set(handles.edit_rename4,'String',handles.sylvalue)

handles=refresh_hist(handles);


guidata(hObject, handles);




function handles=refresh_hist(handles)    

axes(handles.axes_hist);
cla
hist((handles.syltempEnd-handles.syltempStart)/44150,get(handles.slider_classnumber,'Value'))
hold on
plot([get(handles.slider_bound1,'Value') get(handles.slider_bound1,'Value')],[0 10],'red')
plot([get(handles.slider_bound2,'Value') get(handles.slider_bound2,'Value')],[0 10],'green')
plot([get(handles.slider_bound3,'Value') get(handles.slider_bound3,'Value')],[0 10],'black')
hold off
    



function edit_rename1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rename1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rename1 as text
%        str2double(get(hObject,'String')) returns contents of edit_rename1 as a double
duration=(handles.syltempEnd-handles.syltempStart)/44150;
handles.syltempType(duration<=handles.limits(1))=str2double(get(handles.edit_rename1,'String'));
handles.sylType(handles.syltempIndx(duration<=handles.limits(1)))=str2double(get(handles.edit_rename1,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_rename1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rename1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_rename2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rename2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rename2 as text
%        str2double(get(hObject,'String')) returns contents of edit_rename2 as a double
duration=(handles.syltempEnd-handles.syltempStart)/44150;
handles.syltempType((duration<=handles.limits(2)) & (duration>handles.limits(1)))=str2double(get(handles.edit_rename2,'String'));
handles.sylType(handles.syltempIndx((duration<=handles.limits(2)) & (duration>handles.limits(1))))=str2double(get(handles.edit_rename2,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_rename2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rename2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_rename3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rename3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rename3 as text
%        str2double(get(hObject,'String')) returns contents of edit_rename3 as a double
duration=(handles.syltempEnd-handles.syltempStart)/44150;
handles.syltempType((duration<=handles.limits(3)) & (duration>handles.limits(2)))=str2double(get(handles.edit_rename3,'String'));
handles.sylType(handles.syltempIndx((duration<=handles.limits(3)) & (duration>handles.limits(2))))=str2double(get(handles.edit_rename3,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_rename3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rename3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_rename4_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rename4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rename4 as text
%        str2double(get(hObject,'String')) returns contents of edit_rename4 as a double
duration=(handles.syltempEnd-handles.syltempStart)/44150;
handles.syltempType(duration>handles.limits(3))=str2double(get(handles.edit_rename4,'String'));
handles.sylType(handles.syltempIndx(duration>handles.limits(3)))=str2double(get(handles.edit_rename4,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_rename4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rename4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_split1_Callback(hObject, eventdata, handles)
% hObject    handle to slider_split1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

axes(handles.axes_syl1)
hold on
imagesc(handles.MeanSpec(:,:,1))
axis xy
plot([get(handles.slider_split1,'Value') get(handles.slider_split1,'Value')], [0.1 100], 'red')
hold off


guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_split1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_split1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in push_split1.
function push_split1_Callback(hObject, eventdata, handles)
% hObject    handle to push_split1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value=get(handles.slider_split1,'Value')/size(handles.MeanSpec(:,:,1),2);
margin=100;
duration=(handles.syltempEnd-handles.syltempStart)/44150;
crit_bin=(duration<=handles.limits(1));

    a=handles.syltempIndx(crit_bin);
    
for i=1:length(handles.syltempIndx(crit_bin))
    j=a(i)+i-1;
    duration=(handles.syltempEnd-handles.syltempStart)/44150;
    if j==1
        %Variables temporaires
        itemp=1;
        handles.syltempIndx=vertcat(handles.syltempIndx(1:itemp), handles.syltempIndx(itemp)+1, handles.syltempIndx(itemp+1:end)+1);
        valueaddstart=floor(duration(itemp)*44150*value+margin);
        handles.syltempStart=horzcat(handles.syltempStart(1:itemp), handles.syltempStart(itemp)+valueaddstart, handles.syltempStart(itemp+1:end));
        handles.syltempEnd=horzcat(handles.syltempStart(itemp+1)-floor(2*margin), handles.syltempEnd(itemp:end));
        handles.syltempType=vertcat(str2double(get(handles.edit_split_left1,'String')), str2double(get(handles.edit_split_right1,'String')), handles.syltempType(itemp+1:end));
        handles.syltempIndxKey=vertcat(handles.syltempIndxKey(1:itemp), handles.syltempIndxKey(itemp), handles.syltempIndxKey(itemp+1:end));
        
        %variables globale
        handles.sylStart=horzcat(handles.sylStart(1:j), handles.sylStart(j)+valueaddstart, handles.sylStart(j+1:end));
        handles.sylEnd=horzcat(handles.sylStart(j+1)-floor(2*margin), handles.sylEnd(j:end));
        handles.sylType=vertcat(str2double(get(handles.edit_split_left1,'String')), str2double(get(handles.edit_split_right1,'String')), handles.sylType(j+1:end));
        handles.IndxKey=vertcat(handles.IndxKey(1:j), handles.IndxKey(j), handles.IndxKey(j+1:end));
        
    else
        %Variables temporaires
        itemp=find(handles.syltempIndx==j);
        if itemp==1
            handles.syltempIndx=vertcat(handles.syltempIndx(1:itemp), handles.syltempIndx(itemp)+1, handles.syltempIndx(itemp+1:end)+1);
            valueaddstart=floor(duration(itemp)*44150*value+margin);
            handles.syltempStart=horzcat(handles.syltempStart(1:itemp), handles.syltempStart(itemp)+valueaddstart, handles.syltempStart(itemp+1:end));
            handles.syltempEnd=horzcat(handles.syltempStart(itemp+1)-floor(2*margin), handles.syltempEnd(itemp:end));
            handles.syltempType=vertcat(str2double(get(handles.edit_split_left1,'String')), str2double(get(handles.edit_split_right1,'String')), handles.syltempType(itemp+1:end));
            handles.syltempIndxKey=vertcat(handles.syltempIndxKey(1:itemp), handles.syltempIndxKey(itemp), handles.syltempIndxKey(itemp+1:end));
        else
            handles.syltempIndx=vertcat(handles.syltempIndx(1:itemp), handles.syltempIndx(itemp)+1, handles.syltempIndx(itemp+1:end)+1);
            valueaddstart=floor(duration(itemp)*44150*value+margin);
            handles.syltempStart=horzcat(handles.syltempStart(1:itemp), handles.syltempStart(itemp)+valueaddstart, handles.syltempStart(itemp+1:end));
            handles.syltempEnd=horzcat(handles.syltempEnd(1:itemp-1), handles.syltempStart(itemp+1)-floor(2*margin), handles.syltempEnd(itemp:end));
            handles.syltempType=vertcat(handles.syltempType(1:itemp-1), str2double(get(handles.edit_split_left1,'String')), str2double(get(handles.edit_split_right1,'String')), handles.syltempType(itemp+1:end));
            handles.syltempIndxKey=vertcat(handles.syltempIndxKey(1:itemp), handles.syltempIndxKey(itemp), handles.syltempIndxKey(itemp+1:end));
        end
        %variables globale
        handles.sylStart=horzcat(handles.sylStart(1:j), handles.sylStart(j)+valueaddstart, handles.sylStart(j+1:end));
        handles.sylEnd=horzcat(handles.sylEnd(1:j-1), handles.sylStart(j+1)-floor(2*margin), handles.sylEnd(j:end));
        handles.sylType=vertcat(handles.sylType(1:j-1), str2double(get(handles.edit_split_left1,'String')), str2double(get(handles.edit_split_right1,'String')), handles.sylType(j+1:end));
        handles.IndxKey=vertcat(handles.IndxKey(1:j), handles.IndxKey(j), handles.IndxKey(j+1:end));
                
    end
end





guidata(hObject, handles);


function edit_split_left1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_split_left1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_split_left1 as text
%        str2double(get(hObject,'String')) returns contents of edit_split_left1 as a double


% --- Executes during object creation, after setting all properties.
function edit_split_left1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_split_left1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_split_right1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_split_right1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_split_right1 as text
%        str2double(get(hObject,'String')) returns contents of edit_split_right1 as a double


% --- Executes during object creation, after setting all properties.
function edit_split_right1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_split_right1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_split2_Callback(hObject, eventdata, handles)
% hObject    handle to slider_split2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
axes(handles.axes_syl2)
hold on
imagesc(handles.MeanSpec(:,:,2))
axis xy
plot([get(handles.slider_split2,'Value') get(handles.slider_split2,'Value')], [0.1 100], 'red')
hold off


guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_split2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_split2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in push_split2.
function push_split2_Callback(hObject, eventdata, handles)
% hObject    handle to push_split2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value=get(handles.slider_split2,'Value')/size(handles.MeanSpec(:,:,2),2);
margin=100;
duration=(handles.syltempEnd-handles.syltempStart)/44150;
crit_bin=(duration<=handles.limits(2)) & (duration>handles.limits(1));

    a=handles.syltempIndx(crit_bin);
    
for i=1:length(handles.syltempIndx(crit_bin))
    j=a(i)+i-1;
    duration=(handles.syltempEnd-handles.syltempStart)/44150;
    if j==1
        %Variables temporaires
        itemp=1;
        handles.syltempIndx=vertcat(handles.syltempIndx(1:itemp), handles.syltempIndx(itemp)+1, handles.syltempIndx(itemp+1:end)+1);
        valueaddstart=floor(duration(itemp)*44150*value+margin);
        handles.syltempStart=horzcat(handles.syltempStart(1:itemp), handles.syltempStart(itemp)+valueaddstart, handles.syltempStart(itemp+1:end));
        handles.syltempEnd=horzcat(handles.syltempStart(itemp+1)-floor(2*margin), handles.syltempEnd(itemp:end));
        handles.syltempType=vertcat(str2double(get(handles.edit_split_left2,'String')), str2double(get(handles.edit_split_right2,'String')), handles.syltempType(itemp+1:end));
        handles.syltempIndxKey=vertcat(handles.syltempIndxKey(1:itemp), handles.syltempIndxKey(itemp), handles.syltempIndxKey(itemp+1:end));
        
        %variables globale
        handles.sylStart=horzcat(handles.sylStart(1:j), handles.sylStart(j)+valueaddstart, handles.sylStart(j+1:end));
        handles.sylEnd=horzcat(handles.sylStart(j+1)-floor(2*margin), handles.sylEnd(j:end));
        handles.sylType=vertcat(str2double(get(handles.edit_split_left2,'String')), str2double(get(handles.edit_split_right2,'String')), handles.sylType(j+1:end));
        handles.IndxKey=vertcat(handles.IndxKey(1:j), handles.IndxKey(j), handles.IndxKey(j+1:end));
        
    else
        %Variables temporaires
        itemp=find(handles.syltempIndx==j);
        if itemp==1
            handles.syltempIndx=vertcat(handles.syltempIndx(1:itemp), handles.syltempIndx(itemp)+1, handles.syltempIndx(itemp+1:end)+1);
            valueaddstart=floor(duration(itemp)*44150*value+margin);
            handles.syltempStart=horzcat(handles.syltempStart(1:itemp), handles.syltempStart(itemp)+valueaddstart, handles.syltempStart(itemp+1:end));
            handles.syltempEnd=horzcat(handles.syltempStart(itemp+1)-floor(2*margin), handles.syltempEnd(itemp:end));
            handles.syltempType=vertcat(str2double(get(handles.edit_split_left2,'String')), str2double(get(handles.edit_split_right2,'String')), handles.syltempType(itemp+1:end));
            handles.syltempIndxKey=vertcat(handles.syltempIndxKey(1:itemp), handles.syltempIndxKey(itemp), handles.syltempIndxKey(itemp+1:end));
        else
            handles.syltempIndx=vertcat(handles.syltempIndx(1:itemp), handles.syltempIndx(itemp)+1, handles.syltempIndx(itemp+1:end)+1);
            valueaddstart=floor(duration(itemp)*44150*value+margin);
            handles.syltempStart=horzcat(handles.syltempStart(1:itemp), handles.syltempStart(itemp)+valueaddstart, handles.syltempStart(itemp+1:end));
            handles.syltempEnd=horzcat(handles.syltempEnd(1:itemp-1), handles.syltempStart(itemp+1)-floor(2*margin), handles.syltempEnd(itemp:end));
            handles.syltempType=vertcat(handles.syltempType(1:itemp-1), str2double(get(handles.edit_split_left2,'String')), str2double(get(handles.edit_split_right2,'String')), handles.syltempType(itemp+1:end));
            handles.syltempIndxKey=vertcat(handles.syltempIndxKey(1:itemp), handles.syltempIndxKey(itemp), handles.syltempIndxKey(itemp+1:end));
        end
        %variables globale
        handles.sylStart=horzcat(handles.sylStart(1:j), handles.sylStart(j)+valueaddstart, handles.sylStart(j+1:end));
        handles.sylEnd=horzcat(handles.sylEnd(1:j-1), handles.sylStart(j+1)-floor(2*margin), handles.sylEnd(j:end));
        handles.sylType=vertcat(handles.sylType(1:j-1), str2double(get(handles.edit_split_left2,'String')), str2double(get(handles.edit_split_right2,'String')), handles.sylType(j+1:end));
        handles.IndxKey=vertcat(handles.IndxKey(1:j), handles.IndxKey(j), handles.IndxKey(j+1:end));
                
    end
end





guidata(hObject, handles);


function edit_split_left2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_split_left2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_split_left2 as text
%        str2double(get(hObject,'String')) returns contents of edit_split_left2 as a double


% --- Executes during object creation, after setting all properties.
function edit_split_left2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_split_left2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_split_right2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_split_right2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_split_right2 as text
%        str2double(get(hObject,'String')) returns contents of edit_split_right2 as a double


% --- Executes during object creation, after setting all properties.
function edit_split_right2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_split_right2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_split3_Callback(hObject, eventdata, handles)
% hObject    handle to slider_split3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
axes(handles.axes_syl3)
hold on
imagesc(handles.MeanSpec(:,:,3))
axis xy
plot([get(handles.slider_split3,'Value') get(handles.slider_split3,'Value')], [0.1 100], 'red')
hold off


guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_split3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_split3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in push_split3.
function push_split3_Callback(hObject, eventdata, handles)
% hObject    handle to push_split3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value=get(handles.slider_split3,'Value')/size(handles.MeanSpec(:,:,3),2);
margin=100;
duration=(handles.syltempEnd-handles.syltempStart)/44150;
crit_bin=(duration<=handles.limits(3)) & (duration>handles.limits(2));

    a=handles.syltempIndx(crit_bin);
    
for i=1:length(handles.syltempIndx(crit_bin))
    j=a(i)+i-1;
    duration=(handles.syltempEnd-handles.syltempStart)/44150;
    if j==1
        %Variables temporaires
        itemp=1;
        handles.syltempIndx=vertcat(handles.syltempIndx(1:itemp), handles.syltempIndx(itemp)+1, handles.syltempIndx(itemp+1:end)+1);
        valueaddstart=floor(duration(itemp)*44150*value+margin);
        handles.syltempStart=horzcat(handles.syltempStart(1:itemp), handles.syltempStart(itemp)+valueaddstart, handles.syltempStart(itemp+1:end));
        handles.syltempEnd=horzcat(handles.syltempStart(itemp+1)-floor(2*margin), handles.syltempEnd(itemp:end));
        handles.syltempType=vertcat(str2double(get(handles.edit_split_left3,'String')), str2double(get(handles.edit_split_right3,'String')), handles.syltempType(itemp+1:end));
        handles.syltempIndxKey=vertcat(handles.syltempIndxKey(1:itemp), handles.syltempIndxKey(itemp), handles.syltempIndxKey(itemp+1:end));
        
        %variables globale
        handles.sylStart=horzcat(handles.sylStart(1:j), handles.sylStart(j)+valueaddstart, handles.sylStart(j+1:end));
        handles.sylEnd=horzcat(handles.sylStart(j+1)-floor(2*margin), handles.sylEnd(j:end));
        handles.sylType=vertcat(str2double(get(handles.edit_split_left3,'String')), str2double(get(handles.edit_split_right3,'String')), handles.sylType(j+1:end));
        handles.IndxKey=vertcat(handles.IndxKey(1:j), handles.IndxKey(j), handles.IndxKey(j+1:end));
        
    else
        %Variables temporaires
        itemp=find(handles.syltempIndx==j);
        if itemp==1
            handles.syltempIndx=vertcat(handles.syltempIndx(1:itemp), handles.syltempIndx(itemp)+1, handles.syltempIndx(itemp+1:end)+1);
            valueaddstart=floor(duration(itemp)*44150*value+margin);
            handles.syltempStart=horzcat(handles.syltempStart(1:itemp), handles.syltempStart(itemp)+valueaddstart, handles.syltempStart(itemp+1:end));
            handles.syltempEnd=horzcat(handles.syltempStart(itemp+1)-floor(2*margin), handles.syltempEnd(itemp:end));
            handles.syltempType=vertcat(str2double(get(handles.edit_split_left3,'String')), str2double(get(handles.edit_split_right3,'String')), handles.syltempType(itemp+1:end));
            handles.syltempIndxKey=vertcat(handles.syltempIndxKey(1:itemp), handles.syltempIndxKey(itemp), handles.syltempIndxKey(itemp+1:end));
        else
            handles.syltempIndx=vertcat(handles.syltempIndx(1:itemp), handles.syltempIndx(itemp)+1, handles.syltempIndx(itemp+1:end)+1);
            valueaddstart=floor(duration(itemp)*44150*value+margin);
            handles.syltempStart=horzcat(handles.syltempStart(1:itemp), handles.syltempStart(itemp)+valueaddstart, handles.syltempStart(itemp+1:end));
            handles.syltempEnd=horzcat(handles.syltempEnd(1:itemp-1), handles.syltempStart(itemp+1)-floor(2*margin), handles.syltempEnd(itemp:end));
            handles.syltempType=vertcat(handles.syltempType(1:itemp-1), str2double(get(handles.edit_split_left3,'String')), str2double(get(handles.edit_split_right3,'String')), handles.syltempType(itemp+1:end));
            handles.syltempIndxKey=vertcat(handles.syltempIndxKey(1:itemp), handles.syltempIndxKey(itemp), handles.syltempIndxKey(itemp+1:end));
        end
        %variables globale
        handles.sylStart=horzcat(handles.sylStart(1:j), handles.sylStart(j)+valueaddstart, handles.sylStart(j+1:end));
        handles.sylEnd=horzcat(handles.sylEnd(1:j-1), handles.sylStart(j+1)-floor(2*margin), handles.sylEnd(j:end));
        handles.sylType=vertcat(handles.sylType(1:j-1), str2double(get(handles.edit_split_left3,'String')), str2double(get(handles.edit_split_right3,'String')), handles.sylType(j+1:end));
        handles.IndxKey=vertcat(handles.IndxKey(1:j), handles.IndxKey(j), handles.IndxKey(j+1:end));
                
    end
end





guidata(hObject, handles);


function edit_split_left3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_split_left3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_split_left3 as text
%        str2double(get(hObject,'String')) returns contents of edit_split_left3 as a double


% --- Executes during object creation, after setting all properties.
function edit_split_left3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_split_left3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_split_right3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_split_right3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_split_right3 as text
%        str2double(get(hObject,'String')) returns contents of edit_split_right3 as a double


% --- Executes during object creation, after setting all properties.
function edit_split_right3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_split_right3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_split4_Callback(hObject, eventdata, handles)
% hObject    handle to slider_split4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
axes(handles.axes_syl4)
hold on
imagesc(handles.MeanSpec(:,:,4))
axis xy
plot([get(handles.slider_split4,'Value') get(handles.slider_split4,'Value')], [0.1 100], 'red')
hold off


guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_split4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_split4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in push_split4.
function push_split4_Callback(hObject, eventdata, handles)
% hObject    handle to push_split4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value=get(handles.slider_split4,'Value')/size(handles.MeanSpec(:,:,4),2);
margin=100;
duration=(handles.syltempEnd-handles.syltempStart)/44150;
crit_bin=duration>handles.limits(3);

    a=handles.syltempIndx(crit_bin);
    
for i=1:length(handles.syltempIndx(crit_bin))
    j=a(i)+i-1;
    duration=(handles.syltempEnd-handles.syltempStart)/44150;
    if j==1
        %Variables temporaires
        itemp=1;
        handles.syltempIndx=vertcat(handles.syltempIndx(1:itemp), handles.syltempIndx(itemp)+1, handles.syltempIndx(itemp+1:end)+1);
        valueaddstart=floor(duration(itemp)*44150*value+margin);
        handles.syltempStart=horzcat(handles.syltempStart(1:itemp), handles.syltempStart(itemp)+valueaddstart, handles.syltempStart(itemp+1:end));
        handles.syltempEnd=horzcat(handles.syltempStart(itemp+1)-floor(2*margin), handles.syltempEnd(itemp:end));
        handles.syltempType=vertcat(str2double(get(handles.edit_split_left4,'String')), str2double(get(handles.edit_split_right4,'String')), handles.syltempType(itemp+1:end));
        handles.syltempIndxKey=vertcat(handles.syltempIndxKey(1:itemp), handles.syltempIndxKey(itemp), handles.syltempIndxKey(itemp+1:end));
        
        %variables globale
        handles.sylStart=horzcat(handles.sylStart(1:j), handles.sylStart(j)+valueaddstart, handles.sylStart(j+1:end));
        handles.sylEnd=horzcat(handles.sylStart(j+1)-floor(2*margin), handles.sylEnd(j:end));
        handles.sylType=vertcat(str2double(get(handles.edit_split_left4,'String')), str2double(get(handles.edit_split_right4,'String')), handles.sylType(j+1:end));
        handles.IndxKey=vertcat(handles.IndxKey(1:j), handles.IndxKey(j), handles.IndxKey(j+1:end));
        
    else
        %Variables temporaires
        itemp=find(handles.syltempIndx==j);
        if itemp==1
            handles.syltempIndx=vertcat(handles.syltempIndx(1:itemp), handles.syltempIndx(itemp)+1, handles.syltempIndx(itemp+1:end)+1);
            valueaddstart=floor(duration(itemp)*44150*value+margin);
            handles.syltempStart=horzcat(handles.syltempStart(1:itemp), handles.syltempStart(itemp)+valueaddstart, handles.syltempStart(itemp+1:end));
            handles.syltempEnd=horzcat(handles.syltempStart(itemp+1)-floor(2*margin), handles.syltempEnd(itemp:end));
            handles.syltempType=vertcat(str2double(get(handles.edit_split_left4,'String')), str2double(get(handles.edit_split_right4,'String')), handles.syltempType(itemp+1:end));
            handles.syltempIndxKey=vertcat(handles.syltempIndxKey(1:itemp), handles.syltempIndxKey(itemp), handles.syltempIndxKey(itemp+1:end));
        else
            handles.syltempIndx=vertcat(handles.syltempIndx(1:itemp), handles.syltempIndx(itemp)+1, handles.syltempIndx(itemp+1:end)+1);
            valueaddstart=floor(duration(itemp)*44150*value+margin);
            handles.syltempStart=horzcat(handles.syltempStart(1:itemp), handles.syltempStart(itemp)+valueaddstart, handles.syltempStart(itemp+1:end));
            handles.syltempEnd=horzcat(handles.syltempEnd(1:itemp-1), handles.syltempStart(itemp+1)-floor(2*margin), handles.syltempEnd(itemp:end));
            handles.syltempType=vertcat(handles.syltempType(1:itemp-1), str2double(get(handles.edit_split_left4,'String')), str2double(get(handles.edit_split_right4,'String')), handles.syltempType(itemp+1:end));
            handles.syltempIndxKey=vertcat(handles.syltempIndxKey(1:itemp), handles.syltempIndxKey(itemp), handles.syltempIndxKey(itemp+1:end));
        end
        %variables globale
        handles.sylStart=horzcat(handles.sylStart(1:j), handles.sylStart(j)+valueaddstart, handles.sylStart(j+1:end));
        handles.sylEnd=horzcat(handles.sylEnd(1:j-1), handles.sylStart(j+1)-floor(2*margin), handles.sylEnd(j:end));
        handles.sylType=vertcat(handles.sylType(1:j-1), str2double(get(handles.edit_split_left4,'String')), str2double(get(handles.edit_split_right4,'String')), handles.sylType(j+1:end));
        handles.IndxKey=vertcat(handles.IndxKey(1:j), handles.IndxKey(j), handles.IndxKey(j+1:end));
                
    end
end





guidata(hObject, handles);


function edit_split_left4_Callback(hObject, eventdata, handles)
% hObject    handle to edit_split_left4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_split_left4 as text
%        str2double(get(hObject,'String')) returns contents of edit_split_left4 as a double


% --- Executes during object creation, after setting all properties.
function edit_split_left4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_split_left4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_split_right4_Callback(hObject, eventdata, handles)
% hObject    handle to edit_split_right4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_split_right4 as text
%        str2double(get(hObject,'String')) returns contents of edit_split_right4 as a double


% --- Executes during object creation, after setting all properties.
function edit_split_right4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_split_right4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end
