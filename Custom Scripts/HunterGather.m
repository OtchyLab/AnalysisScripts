function varargout = HunterGather(varargin)
% HUNTERGATHER M-file for HunterGather.fig
%      HUNTERGATHER, by itself, creates a new HUNTERGATHER or raises the existing
%      singleton*.
%
%      H = HUNTERGATHER returns the handle to a new HUNTERGATHER or the handle to
%      the existing singleton*.
%
%      HUNTERGATHER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HUNTERGATHER.M with the given input arguments.
%
%      HUNTERGATHER('Property','Value',...) creates a new HUNTERGATHER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HunterGather_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HunterGather_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HunterGather

% Last Modified by GUIDE v2.5 30-Jun-2010 18:47:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HunterGather_OpeningFcn, ...
                   'gui_OutputFcn',  @HunterGather_OutputFcn, ...
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

% --- Executes just before HunterGather is made visible.
function HunterGather_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HunterGather (see VARARGIN)

% Choose default command line output for HunterGather
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes HunterGather wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = HunterGather_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


function edit_DOB_Callback(hObject, eventdata, handles)
% hObject    handle to edit_DOB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_DOB as text
%        str2double(get(hObject,'String')) returns contents of edit_DOB as a double

handles.DOB = get(handles.edit_DOB, 'String');

splits = regexp(handles.DOB, '-', 'split'); %Parse file name for date
year = str2double(char(splits{1}(3)));
month = str2double(char(splits{1}(1)));
day = str2double(char(splits{1}(2)));

handles.DOB = datenum(year, month, day);

guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function edit_DOB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_DOB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in push_AddSim.
function push_AddSim_Callback(hObject, eventdata, handles)
% hObject    handle to push_AddSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Request user input for where to find next Similarity File
[handles.SimFileName, handles.SimFileLoc]=uigetfile('*.mat','Select Sim File:');
splits = regexp(handles.SimFileName, '_', 'split'); %Parse file name for date
year = str2double(['20' char(splits(2))]);
month = str2double(char(splits(3)));
day = str2double(char(splits{4}(1:2)));
stringDate = [['20' char(splits(2))] '-' char(splits(3)) '-' char(splits{4}(1:2))];

date = datenum(year, month, day);
DPH = date - handles.DOB;

DPHstring = num2str(DPH);

load([handles.SimFileLoc '\' handles.SimFileName]);
if ~isfield(handles,'Var')
    handles.Var = [];
    handles.VarDPH = [];
    handles.SimList = [];
end
handles.Var = [handles.Var; SimAnal.SetVar];
handles.VarDPH = [handles.VarDPH; DPH];
handles.SimList = [handles.SimList, stringDate '  -  ' DPHstring 'DPH' sprintf('\n')];

axes(handles.axes_List);
cla;
text(0, 1, handles.SimList);
set(handles.axes_List,'Visible','off');
box on

guidata(hObject, handles);

end


% --- Executes on button press in push_PlotVar.
function push_PlotVar_Callback(hObject, eventdata, handles)
% hObject    handle to push_PlotVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot(handles.VarDPH,handles.Var);
xlim([35 max(SylSimAll(1,1))+15]);
ylim([0 100]);
ylabel('Variability (ADU)');
xlabel('Days Post Hatch');
hold on;

end

% --- Executes on button press in push_SetDirectory.
function push_SetDirectory_Callback(hObject, eventdata, handles)
% hObject    handle to push_SetDirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[folder]=uigetdir('Select Parent Folder:');
cd(folder);

guidata(hObject, handles);

end

% --- Executes on button press in push_ExportVar.
function push_ExportVar_Callback(hObject, eventdata, handles)
% hObject    handle to push_ExportVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end
