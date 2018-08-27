function varargout = tiffAngleHelper(varargin)
% TIFFANGLEHELPER MATLAB code for tiffAngleHelper.fig
%      TIFFANGLEHELPER, by itself, creates a new TIFFANGLEHELPER or raises the existing
%      singleton*.
%
%      H = TIFFANGLEHELPER returns the handle to a new TIFFANGLEHELPER or the handle to
%      the existing singleton*.
%
%      TIFFANGLEHELPER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TIFFANGLEHELPER.M with the given input arguments.
%
%      TIFFANGLEHELPER('Property','Value',...) creates a new TIFFANGLEHELPER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tiffAngleHelper_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tiffAngleHelper_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tiffAngleHelper

% Last Modified by GUIDE v2.5 25-Jul-2018 15:17:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tiffAngleHelper_OpeningFcn, ...
                   'gui_OutputFcn',  @tiffAngleHelper_OutputFcn, ...
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


% --- Executes just before tiffAngleHelper is made visible.
function tiffAngleHelper_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tiffAngleHelper (see VARARGIN)

% Do all the stuff
% select file manually

%%%%%%%%%%%%%%%%%%%%%%%%%
    file = 'TestCube_V-1.06_00001_00001.tif';
    filepath = strcat('/Users/sadiela/Desktop/TestCubeTiffs/', file);
    global tiff;
 
    tiff = [];
    i = 0;
    %title(handles.axes1, 'Loading image...');
    drawnow;
    try
        while true
            i = i + 1;
            tiff(i,:,:) = imread(filepath, i); % dis might be d issue: handles.files is just a cell array of strings
        end
    catch ME
    end
            
    % set slider values for this file
    [handles.numIm, ~, ~] = size(tiff);
    set(handles.zslider, "Min", 1);
    set(handles.zslider, "Max", handles.numIm);
    set(handles.zslider, "SliderStep", [1/(handles.numIm - 1), 0.10]);
    handles.edit_stepSize.String = "1";
    set(handles.zslider, "Value", 1);


    zslider_Callback(handles.zslider, [], handles);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose default command line output for tiffAngleHelper
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tiffAngleHelper wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = tiffAngleHelper_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function [x, y] = getScatter(handles)
    binIm = handles.binIm;
    [rows, cols] = size(binIm);
    points = zeros(2, 1); % this will hold the x&y coordinates for the first "white' pixel in each row
    counter = 1;
    for i = 1:rows
        col = 1;
        while (binIm(i, col) == 1 && col < cols)
            col = col+1;
        end
        if (col ~= cols) 
            points(1,counter) = i; % y components (row #)
            points(2,counter) = col; % x components (column #)
            counter = counter+1;
        %else
            %points(1,i) = 0;
            %points(2,i) = 0;
        end
    end

    x = fliplr(points(2,:));
    y = points(1,:);
    % get the scatter data for linreg
    % step 1: get current binary data


   
function edit_goto_Callback(hObject, eventdata, handles)

function edit_goto_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_stepSize_Callback(hObject, eventdata, handles)
size = str2num(handles.edit_stepSize.String);
set(handles.zslider, "SliderStep", [size/(handles.numIm - 1), size/(handles.numIm - 1)]);
guidata(hObject, handles);

function edit_stepSize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function zslider_Callback(hObject, eventdata, handles)
 global tiff;
    
    ulim = size(tiff, 1);
    colormap bone;
    v = get(handles.zslider, 'Value');
    %imNumSelected = int32(get(handles.zslider, "Value"));
    
    imagesc(squeeze(tiff(v,:,:)), 'Parent', handles.axes1);
    title(handles.axes1, sprintf('Slice %d of %d', v, ulim));
    drawnow;

function zslider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in push_go.
function push_go_Callback(hObject, eventdata, handles)
handles.pos = str2num(handles.edit_goto.String);
handles.zslider.Value = handles.pos;
guidata(hObject, handles);
zslider_Callback(handles.zslider, [], handles);


% --- Executes on button press in push_binary.
function push_binary_Callback(hObject, eventdata, handles)
global tiff;
% get binary image data for the plot
% step 1: get current tiff
v = handles.zslider.Value;
binTiff = squeeze(tiff(v,11:502,11:502));
for j = 1:492
    for k = 1:492
        if(binTiff(j,k) <= 5500)
            binTiff(j,k) = 0;
        else
            binTiff(j,k) = 1;
        end
    end
end
handles.binIm = binTiff;
guidata(hObject, handles);
set(handles.axes2, 'Units', 'pixels');
resizePos = get(handles.axes2, 'Position');
binTiff = imresize(binTiff, [resizePos(3) resizePos(3)]);
axes(handles.axes2);
imshow(binTiff);
set(handles.axes2, 'Units', 'normalized');
guidata(hObject, handles);



% --- Executes on button press in push_linear.
function push_linear_Callback(hObject, eventdata, handles)
axes(handles.axes3);
[x, y] = getScatter(handles);
[xreg, yreg] = simpleLinReg(x, y);
scatter(x, y);
hold;
plot(xreg,yreg);
handles.axes3.XLim = [0, 500];
handles.axes3.YLim = [0, 500];
