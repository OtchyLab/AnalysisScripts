function varargout = Display_spectrogram(varargin)
% DISPLAY_SPECTROGRAM M-file for Display_spectrogram.fig
%      DISPLAY_SPECTROGRAM, by itself, creates a new DISPLAY_SPECTROGRAM or raises the existing
%      singleton*.
%
%      H = DISPLAY_SPECTROGRAM returns the handle to a new DISPLAY_SPECTROGRAM or the handle to
%      the existing singleton*.
%
%      DISPLAY_SPECTROGRAM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISPLAY_SPECTROGRAM.M with the given input arguments.
%
%      DISPLAY_SPECTROGRAM('Property','Value',...) creates a new DISPLAY_SPECTROGRAM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Display_spectrogram_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Display_spectrogram_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help Display_spectrogram

% Last Modified by GUIDE v2.5 21-Dec-2005 12:18:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Display_spectrogram_OpeningFcn, ...
                   'gui_OutputFcn',  @Display_spectrogram_OutputFcn, ...
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


% --- Executes just before Display_spectrogram is made visible.
function Display_spectrogram_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Display_spectrogram (see VARARGIN)
set(hObject, 'Units', 'pixels');



% Choose default command line output for Display_spectrogram
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Display_spectrogram wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Display_spectrogram_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function Spectrogram_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Spectrogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate Spectrogram


% --- Executes on button press in SyllA.
function SyllA_Callback(hObject, eventdata, handles)
% hObject    handle to SyllA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in SyllB.
function SyllB_Callback(hObject, eventdata, handles)
% hObject    handle to SyllB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in SyllC.
function SyllC_Callback(hObject, eventdata, handles)
% hObject    handle to SyllC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in SyllD.
function SyllD_Callback(hObject, eventdata, handles)
% hObject    handle to SyllD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in SyllE.
function SyllE_Callback(hObject, eventdata, handles)
% hObject    handle to SyllE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in SyllF.
function SyllF_Callback(hObject, eventdata, handles)
% hObject    handle to SyllF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


