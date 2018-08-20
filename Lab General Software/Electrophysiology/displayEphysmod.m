function varargout = displayEphys(varargin)
% displayEphys M-file for displayEphy.fig
%      displayEphys, by itself, creates a new displayEphys or raises the existing
%      singleton*.
%
%      H = displayEphys returns the handle to a new displayEphys or the handle to
%      the existing singleton*.
%
%      displayEphys('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in displayEphys.M with the given input arguments.
%
%      displayEphys('Property','Value',...) creates a new displayEphys or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before test_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to test_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help test

% Last Modified by GUIDE v2.5 20-Jun-2006 20:12:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @test_OpeningFcn, ...
                   'gui_OutputFcn',  @test_OutputFcn, ...
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


% --- Executes just before test is made visible.
function test_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to test (see VARARGIN)

% Choose default command line output for test
handles.output = hObject;
% Update handles structure
handles.channel=1;
handles.counter=1;
handles.timestamp=1;
handles.motif_counter=1;
guidata(hObject, handles);
linkaxes([handles.spectrogram handles.audio handles.electrode handles.rasters],'x');
xlabel(handles.electrode,'Time(s)');

% UIWAIT makes test wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = test_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in channel.
function channel_Callback(hObject, eventdata, handles)
% hObject    handle to channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns channel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from channel
handles.channel=get(hObject,'Value');
DisplaySignal(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.filename,handles.dirpath]=uigetfile('*.*');
load([handles.dirpath handles.filename])
cd (handles.dirpath);
exper.dir=handles.dirpath;
handles.ch=exper.sigCh;
handles.fs=exper.desiredInSampRate;
handles.audioCh=exper.audioCh;
handles.exper=exper;
handles.counter=1;
handles.ncounter=1;

guidata(hObject, handles);


function recnum_Callback(hObject, eventdata, handles)
% hObject    handle to recnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recnum as text
handles.rec=str2double(get(hObject,'String'));
set(handles.recnum,'String',num2str(handles.rec));
DisplayData(hObject,eventdata,handles);
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function recnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% function handlesOut = DisplayData(handles)
function DisplayData(hObject,eventData,handles)

%display data

[handles.amplitude, time, HWChannels, startSamp, timeCreated] = loadData(handles.exper,handles.rec,handles.audioCh);
[handles.signal, time, HWChannels, startSamp, timeCreated] = loadData(handles.exper,handles.rec,handles.ch(handles.channel));
subplot(handles.spectrogram); 
specgram1(handles.amplitude,512,handles.fs,400,350);ylim([0 10000]);
subplot(handles.audio);
plot(0:1/handles.fs:(length(handles.amplitude)/handles.fs)-1/handles.fs,handles.amplitude);xlim([0 (length(handles.amplitude)/handles.fs)-1/handles.fs]);ylim([min(handles.amplitude) max(handles.amplitude)]);
subplot(handles.electrode);
plot(0:1/handles.fs:(length(handles.signal)/handles.fs)-1/handles.fs,handles.signal);xlim([0 (length(handles.signal)/handles.fs)-1/handles.fs]);ylim([min(handles.signal) max(handles.signal)]);
subplot(handles.rasters);
plot(0:0.01:0.02,0:0.01:0.02);xlim([0 (length(handles.signal)/handles.fs)-1/handles.fs]);ylim([-1 4]);
xlabel(handles.rasters,'Time(s)');
handlesOut=handles;
[handles.syllStartTimes, handles.syllEndTimes, noiseEst, noiseStd, soundEst, soundStd] = aSAP_segSyllablesFromRawAudio(handles.amplitude, handles.fs);
subplot(handles.spectrogram);
for i=1:length(handles.syllStartTimes)
   line([handles.syllStartTimes(i) handles.syllStartTimes(i)],[0 10000],'color','w');
end
for i=1:length(handles.syllStartTimes)
   line([handles.syllEndTimes(i) handles.syllEndTimes(i)],[0 10000],'color','y');
end
% handles.spectimg = findobj('Type','image');
% set(handles.spectimg,'ButtonDownFcn',{mousepress(handles.counter,handles.timestamp,handles.ncounter)});
guidata(hObject, handles);


function DisplaySignal(handles)
[handles.signal, time, HWChannels, startSamp, timeCreated] = loadData(handles.exper,handles.rec,handles.ch(handles.channel));
subplot(handles.electrode);
plot(0:1/handles.fs:(length(handles.signal)/handles.fs)-1/handles.fs,handles.signal);xlim([0 (length(handles.signal)/handles.fs)-1/handles.fs]);ylim([min(handles.signal) max(handles.signal)]);
xlabel(handles.electrode,'Time(s)');


% --- Executes on button press in Motif.
function Motif_Callback(hObject, eventdata, handles)
% hObject    handle to Motif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.counter=1;
    handles.motif_counter=handles.motif_counter+1;
    htext=text(mean(handles.timestamp),8000,['Motif ' num2str(handles.motif_counter)],'HorizontalAlignment','center','Color',[1 1 0]);
 

guidata(hObject, handles);    



% function mousepress(counter,timestamp,ncounter)   
%     if (counter<3)    
%         click = get(gca,'CurrentPoint')
%         timestamp(counter)=click(1,1);
%         line([timestamp(counter) timestamp(counter)],[0 10000],'Color','r');
%        
%     else
%         ncounter=abs(mod(counter,2)-2);
%         click = get(gca,'CurrentPoint')
%         timestamp(ncounter)=click(1,1);
%         line([timestamp(ncounter) timestamp(ncounter)],[0 10000], 'Color','r');
%     end
% counter=counter+1;

