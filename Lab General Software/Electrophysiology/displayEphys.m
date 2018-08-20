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

% Last Modified by GUIDE v2.5 03-Sep-2006 11:15:10
% This GUI displays sound and ephys data and allows
% the user to define motifs. 

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
%Initialize the variables
% handles.threshold=1; %default thresholf for syllable segmentation
% handles.xlimit=0;
% handles.motif_no=0;
% handles.save=0;
% handles.counter=1;
% handles.spikemarker=0;
% set(handles.thresh,'String',num2str(handles.threshold));
% handles.channel=1;
%Link axes together
linkaxes([handles.spectrogram handles.audio handles.electrode handles.rasters],'x');
xlabel(handles.electrode,'Time(s)');
%define cell 
handles.cell.ID='cell';
handles.cell.rec=0;
% handles.cell.motif=[0 0];
% handles.cell.syll=[0 0];
% handles.cell.spikes=[0 0];
guidata(hObject, handles);
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

%Determines which channel(electrode to display)
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

%Load in the experiment 
[handles.filename,handles.dirpath]=uigetfile('*.*');
load([handles.dirpath handles.filename]);
cd (handles.dirpath);
exper.dir=handles.dirpath;
handles.ch=exper.sigCh;
if (exper.datecreated=='0')
    handles.alex=1; % define alex's data this way. also look at LoadDataAlex
else
    handles.alex=0;
end
handles.spikes=0;
handles.fs=exper.desiredInSampRate;
handles.audioCh=exper.audioCh;
handles.exper=exper;
handles.birdname=exper.birdname;
%Initialize the parameters
handles.threshold=1; %default thresholf for syllable segmentation
handles.xlimit=0;
handles.motif_no=0;
handles.song_no=0;
handles.save=0;
handles.counter=1;
handles.spikemarker=0;
set(handles.thresh,'String',num2str(handles.threshold));
handles.channel=1;
guidata(hObject, handles);



function recnum_Callback(hObject, eventdata, handles)
% hObject    handle to recnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recnum as text
handles.rec=str2double(get(hObject,'String'));
set(handles.recnum,'String',num2str(handles.rec));
handles.xlimit=0;
handles=DisplayData(handles);
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

function handlesOut = DisplayData(handles)
%display data
handles.spiketime=0;
handles.counter=1;
% handles.spikes=0;
%alex's data format
if (handles.alex==1)
    handles.amplitude=loadDataAlex(handles.rec, handles.audioCh);
    handles.signal=loadDataAlex(handles.rec, handles.ch(handles.channel));
else
%aaron's data format
 [handles.amplitude, time, HWChannels, startSamp, timeCreated] = loadData(handles.exper,handles.rec,handles.audioCh);
 [handles.signal, time, HWChannels, startSamp, timeCreated] = loadData(handles.exper,handles.rec,handles.ch(handles.channel));
end
subplot(handles.spectrogram); 
cla;
%specgram_hc(handles.amplitude,512,handles.fs,400,350);ylim([0 10000]);
 specgram1(handles.amplitude,512,handles.fs,400,350);ylim([0 10000]);
%Find the hand;le for the spectrogram
specthandle = findobj('Type','image');
%Define the function to execute whenever the mouse clicks on the
%spectrogram
set(specthandle,'ButtonDownFcn','displayEphys(''spectrogram_Callback'',gcbo,[],guidata(gcbo))');
subplot(handles.audio);
cla;
plot(0:1/handles.fs:(length(handles.amplitude)/handles.fs)-1/handles.fs,handles.amplitude);ylim([min(handles.amplitude) max(handles.amplitude)]);
subplot(handles.electrode);
cla;
plot(0:1/handles.fs:(length(handles.signal)/handles.fs)-1/handles.fs,handles.signal);
[handles.syllStartTimes, handles.syllEndTimes, noiseEst, noiseStd, soundEst, soundStd] = aSAP_segSyllablesFromRawAudio(handles.amplitude, handles.fs, handles.threshold);
subplot(handles.spectrogram);
for i=1:length(handles.syllStartTimes)
   handles.linestart(i)=line([handles.syllStartTimes(i) handles.syllStartTimes(i)],[0 10000],'color','w');
end
for i=1:length(handles.syllStartTimes)
   handles.lineend(i)=line([handles.syllEndTimes(i) handles.syllEndTimes(i)],[0 10000],'color','y');
end
if(handles.xlimit==0)
    xlimitstart=max(handles.syllStartTimes(1)-0.3,0);
    xlimitend=min(length(handles.amplitude)/handles.fs,(handles.syllEndTimes(end)+0.3));
    handles.xlimit=([xlimitstart xlimitend]);
end
xlim(handles.xlimit)

subplot(handles.rasters);
cla;
%load spikes
if(handles.spikes==1)
    filename = getspikeDatafile(handles.exper, handles.rec,handles.ch(handles.channel));
    load([filename]);
    cluster=cluster_class;
    cluster(:,2)=cluster(:,2)/1000;
    cellindex=find(cluster(:,1)==1);
    handles.sorted_spikes=cluster(cellindex,2);
    cellindex=find(cluster(:,1)==0);
    handles.unsorted_spikes=cluster(cellindex,2);
    subplot(handles.rasters);
% axis off;
        for i=1:length(handles.sorted_spikes)
            line([handles.sorted_spikes(i) handles.sorted_spikes(i)],[0.3 1],'color','r');
        end
    hold on;
    y = zeros(length(handles.unsorted_spikes));
    plot(handles.unsorted_spikes,y,'b+');ylim([-0.3 1.3]);
    hold off;
end
handlesOut=handles;


function DisplaySignal(handles)
    
if (handles.alex==1)
    handles.signal=loadDataAlex(handles.rec, handles.ch(handles.channel));
else
    [handles.signal, time, HWChannels, startSamp, timeCreated] = loadData(handles.exper,handles.rec,handles.ch(handles.channel));
end
subplot(handles.electrode);
handles.xlimit=get(gca,'xlim');
plot(0:1/handles.fs:(length(handles.signal)/handles.fs)-1/handles.fs,handles.signal);ylim([min(handles.signal) max(handles.signal)]);
xlabel(handles.electrode,'Time(s)');
xlim(handles.xlimit);



% --- Executes during object creation, after setting all properties.
function spectrogram_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spectrogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate spectrogram

function spectrogram_Callback(hObject, eventdata, handles)
% hObject    handle to recnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.counter<3)   
    click = get(gca,'CurrentPoint');
    handles.markertime(handles.counter)=click(1,1);
    handles.marker(handles.counter) = line([handles.markertime(handles.counter) handles.markertime(handles.counter)],[0 10000],'color','r'); 
else
    ncounter=abs(rem(handles.counter,2)-2);
    delete(handles.marker(ncounter)); 
    click = get(gca,'CurrentPoint');
    handles.markertime(ncounter)=click(1,1);
    handles.marker(ncounter) = line([handles.markertime(ncounter) handles.markertime(ncounter)],[0 10000],'colo','r');
end
handles.counter=handles.counter+1;   
guidata(hObject, handles);

function rasters_Callback(hObject, eventdata, handles)
% hObject    handle to recnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
click = get(gca,'CurrentPoint');
if (handles.spiketime~=0)
    delete(handles.spikemarker);
end
handles.spiketime=click(1,1);
hold on
handles.spikemarker = plot(handles.spiketime, -0.1,'g:*'); 
hold off
guidata(hObject, handles);

% --- Executes on button press in previous.
function previous_Callback(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.rec=handles.rec-1;
handles.xlimit=0;
handles=DisplayData(handles);
set(handles.recnum,'String',num2str(handles.rec));
guidata(hObject, handles);

% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.rec=handles.rec+1;
handles.xlimit=0;
handles=DisplayData(handles);
set(handles.recnum,'String',num2str(handles.rec));
guidata(hObject, handles);


% --- Executes on button press in clusters.
function clusters_Callback(hObject, eventdata, handles)
% hObject    handle to clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 if (handles.alex==1)
     if (handles.rec<10)
        fill='000';
     elseif (handles.rec<100)
        fill='00';
     else
        fill='0';
     end
     filename=['times_adult_0714D-f' fill num2str(handles.rec) '.mat']; 
 else
    filename = getspikeDatafile(handles.exper, handles.rec,handles.ch(handles.channel));
 end
 load([filename]);
 cluster=cluster_class;
 cluster(:,2)=cluster(:,2)/1000;
 cellindex=find(cluster(:,1)==1);
 handles.sorted_spikes=cluster(cellindex,2);
 cellindex=find(cluster(:,1)==0);
 handles.unsorted_spikes=cluster(cellindex,2);
 subplot(handles.rasters);
 cla;
 handles.xlimit=get(gca,'xlim');
% axis off;
     for i=1:length(handles.sorted_spikes)
         line([handles.sorted_spikes(i) handles.sorted_spikes(i)],[0.3 1],'color','r');
     end
 hold on;
 y = zeros(length(handles.unsorted_spikes));
 plot(handles.unsorted_spikes,y,'b+');xlim(handles.xlimit);ylim([-0.3 1.3]);
 hold off;

handles.spikes=1;
guidata(hObject, handles);



% --- Executes on slider movement.
function thresh_Callback(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.threshold=str2double(get(hObject,'String'));
set(handles.thresh,'String',num2str(handles.threshold));
for i=1:length(handles.syllStartTimes)
    delete(handles.linestart(i));
    delete(handles.lineend(i));
end
[handles.syllStartTimes, handles.syllEndTimes, noiseEst, noiseStd, soundEst, soundStd] = aSAP_segSyllablesFromRawAudio(handles.amplitude, handles.fs, handles.threshold);
subplot(handles.spectrogram);

for i=1:length(handles.syllStartTimes)
    handles.lineend(i)=line([handles.syllEndTimes(i) handles.syllEndTimes(i)],[0 10000],'color','y');
    handles.linestart(i)=line([handles.syllStartTimes(i) handles.syllStartTimes(i)],[0 10000],'color','w');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on button press in delete.
function delete_Callback(hObject, eventdata, handles)
% hObject    handle to delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.counter>1)
    counter=abs(rem(handles.counter-1,2)-2);
    index=find(handles.syllEndTimes > handles.markertime(counter));
    handles.syllEndTimes=[handles.syllEndTimes(1:index(1)-1),handles.syllEndTimes(index(1)+1:end)];
    handles.syllStartTimes=[handles.syllStartTimes(1:index(1)-1),handles.syllStartTimes(index(1)+1:end)];
    delete(handles.linestart(index(1)));delete(handles.lineend(index(1)));
    handles.linestart=[handles.linestart(1:index(1)-1),handles.linestart(index(1)+1:end)];
    handles.lineend=[handles.lineend(1:index(1)-1),handles.lineend(index(1)+1:end)];
end
guidata(hObject, handles);

% --- Executes on button press in startsyll. Add the start of a new syllable
% that was not recognized by the program
function startsyll_Callback(hObject, eventdata, handles)
% hObject    handle to startsyll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%insert a new start time into syllstarttimes, draw a line, add this to line
%start and increment the pointers to the other lines
if (handles.counter>1)
    counter=abs(rem(handles.counter-1,2)-2);
    index=find(handles.syllStartTimes > handles.markertime(counter));
    if (numel(index)==0)
        handles.syllStartTimes=[handles.syllStartTimes,handles.markertime(counter)];
        handles.linestart(end+1)=line([handles.markertime(counter) handles.markertime(counter)],[0 10000],'color','w');
    else
        handles.syllStartTimes=[handles.syllStartTimes(1:index(1)-1),handles.markertime(counter),handles.syllStartTimes(index(1):end)];
        temp_linehandle=handles.linestart(index(1):end);
        handles.linestart(index(1))=line([handles.markertime(counter) handles.markertime(counter)],[0 10000],'color','w');
        handles.linestart=[handles.linestart(1:index(1)),temp_linehandle];
    end
end
 guidata(hObject, handles);
 
% --- Executes on button press in endsyll. Add the end of a new syllable
% that was not recognized by the program
function endsyll_Callback(hObject, eventdata, handles)
% hObject    handle to startsyll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.counter>1)
    counter=abs(rem(handles.counter-1,2)-2);
    index=find(handles.syllEndTimes > handles.markertime(counter));
    if (numel(index)==0)
        handles.syllEndTimes=[handles.syllEndTimes,handles.markertime(counter)];
        handles.lineend(end+1)=line([handles.markertime(counter) handles.markertime(counter)],[0 10000],'color','w');
    else
        handles.syllEndTimes=[handles.syllEndTimes(1:index(1)-1),handles.markertime(counter),handles.syllEndTimes(index(1):end)];
        temp_linehandle=handles.lineend(index(1):end);
        handles.lineend(index(1))=line([handles.markertime(counter) handles.markertime(counter)],[0 10000],'color','y');
        handles.lineend=[handles.lineend(1:index(1)),temp_linehandle];
    end
end
guidata(hObject, handles);


% --- Executes on button press in savecell.
function savecell_Callback(hObject, eventdata, handles)
% hObject    handle to savecell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.motif_no=0;
handles.song_no=0;
handles.save=1;
handles.cell.rec=0;
handles.cell.motif=[0 0];
handles.cell.syll=[0 0];
handles.cell.spikes=[0 0];
[handles.cellName,handles.cellPath]=uiputfile('cell.mat','Save cell as:');
presentCell=handles.cell;
save ([handles.cellPath,handles.cellName], 'presentCell');
guidata(hObject, handles);

% --- Executes on button press in addmotif.
function addmotif_Callback(hObject, eventdata, handles)
% hObject    handle to addmotif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (length(handles.markertime)>1 & handles.save==1 & handles.spikes==1)
    button='Go Ahead';
    handles.motif_no=handles.motif_no+1;
    ncounter=abs(rem(handles.counter-1,2)-2)
    delete(handles.marker(1)); 
    delete(handles.marker(2));
    startMotif=min(handles.markertime(1),handles.markertime(2))
    endMotif=max(handles.markertime(1),handles.markertime(2))
    startMotifIndex=find(handles.syllStartTimes>startMotif,1,'first')
    endMotifIndex=find(handles.syllEndTimes<endMotif,1,'last')
    handles.cell.ID=handles.birdname;
    handles.cell.rec(handles.motif_no)=handles.rec;
    handles.cell.motif(handles.motif_no,1)=handles.syllStartTimes(startMotifIndex);
    handles.cell.motif(handles.motif_no,2)=handles.syllEndTimes(endMotifIndex);
    
    no_syll=find(handles.syllStartTimes>startMotif & handles.syllStartTimes<endMotif);
    no_syll=length(no_syll);
    if (handles.motif_no>1 & no_syll~=handles.no_syll)
        button = questdlg('The number of syllables is different from previous motifs','Warning!','New','Go Ahead','New');
    end
    if (handles.motif_no>1 & strcmp(button, 'New'))
            handles.motif_no=handles.motif_no-1;
    else
    
        handles.no_syll=no_syll;


        for i=1:2:(2*no_syll)
            handles.cell.syll(handles.motif_no,i)=handles.syllStartTimes(startMotifIndex+floor(i/2));
            handles.cell.syll(handles.motif_no,i+1)=handles.syllEndTimes(startMotifIndex+floor(i/2));
        end
        if (handles.spikes==1)
            spike_index=find(handles.sorted_spikes>handles.syllStartTimes(startMotifIndex)-0.3 & handles.sorted_spikes<handles.syllEndTimes(endMotifIndex)+0.3);
            handles.cell.spikes(handles.motif_no,1:(spike_index(end)-spike_index(1))+1)=handles.sorted_spikes(spike_index(1):spike_index(end));
        end
        presentCell=handles.cell;
        save ([handles.cellPath,handles.cellName], 'presentCell');
        %add a screen pointer to the motif that has already been selected
        subplot(handles.spectrogram);
        for i=startMotifIndex:(no_syll-1) + startMotifIndex    
            delete(handles.linestart(i));
            delete(handles.lineend(i));
        end
        handles.motif_line=line([handles.syllStartTimes(startMotifIndex) handles.syllEndTimes(endMotifIndex)],[6000 6000],'color','r');
        handles.motif_text=text((handles.syllStartTimes(startMotifIndex)+handles.syllEndTimes(endMotifIndex))/2,8000,['Motif Number ' num2str(handles.motif_no)],'HorizontalAlignment','center','Color',[1 1 0]);

        if (endMotifIndex==length(handles.syllEndTimes))
            handles.syllEndTimes=handles.syllEndTimes(1:startMotifIndex-1);
            handles.syllStartTimes=handles.syllStartTimes(1:startMotifIndex-1);
            handles.linestart=handles.linestart(1:startMotifIndex-1);
            handles.lineend=handles.lineend(1:startMotifIndex-1);
        else    
            handles.syllEndTimes=[handles.syllEndTimes(1:startMotifIndex-1),handles.syllEndTimes(endMotifIndex+1:end)];
            handles.syllStartTimes=[handles.syllStartTimes(1:startMotifIndex-1),handles.syllStartTimes(endMotifIndex+1:end)];
            handles.linestart=[handles.linestart(1:startMotifIndex-1),handles.linestart(endMotifIndex+1:end)];
            handles.lineend=[handles.lineend(1:startMotifIndex-1),handles.lineend(endMotifIndex+1:end)];
        end
    
    end
    handles.counter=1
    guidata(hObject, handles);
else
    warndlg('Load spikes (or cell) to add Motif')
end

% --- Executes on button press in deletespike.
function addspike_Callback(hObject, eventdata, handles)
% hObject    handle to deletespike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.spiketime~=0)
    added_spike1 = find (handles.unsorted_spikes<handles.spiketime, 1,'last');
    added_spike2 = find (handles.unsorted_spikes>handles.spiketime, 1,'first');
    if (abs(handles.unsorted_spikes(added_spike1)-handles.spiketime)<abs(handles.unsorted_spikes(added_spike2)-handles.spiketime))
        added_spike=added_spike1;
    else
        added_spike=added_spike2;
    end
    added_spike_index=find(handles.sorted_spikes<handles.unsorted_spikes(added_spike), 1,'last');
    handles.sorted_spikes=[handles.sorted_spikes(1:added_spike_index); handles.unsorted_spikes(added_spike); handles.sorted_spikes((added_spike_index+1):end)];
    handles.unsorted_spikes=[handles.unsorted_spikes(1:added_spike-1); handles.unsorted_spikes(added_spike+1:end)];
    handles.spiketime=0;
    subplot(handles.rasters);
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
guidata(hObject, handles);


% --- Executes on button press in deletespike.
function deletespike_Callback(hObject, eventdata, handles)
% hObject    handle to deletespike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.spiketime~=0)
    added_spike1 = find (handles.sorted_spikes<handles.spiketime, 1,'last');
    added_spike2 = find (handles.sorted_spikes>handles.spiketime, 1,'first');
    if (abs(handles.sorted_spikes(added_spike1)-handles.spiketime)<abs(handles.sorted_spikes(added_spike2)-handles.spiketime));
        added_spike=added_spike1;
    else
        added_spike=added_spike2;
    end
    added_spike_index=find(handles.unsorted_spikes<handles.sorted_spikes(added_spike), 1,'last');
    handles.unsorted_spikes=[handles.unsorted_spikes(1:added_spike_index); handles.sorted_spikes(added_spike); handles.unsorted_spikes((added_spike_index+1):end)];
    handles.sorted_spikes=[handles.sorted_spikes(1:added_spike-1); handles.sorted_spikes(added_spike+1:end)];
    handles.spiketime=0;
    
    subplot(handles.rasters);
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
guidata(hObject, handles);

% --- Executes on button press in scanback.
function scanback_Callback(hObject, eventdata, handles)
% hObject    handle to scanback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xlimit=get(gca,'xlim');
handles.xlimit=[2*xlimit(1)-xlimit(2) xlimit(1)];
xlim(handles.xlimit);
guidata(hObject, handles);

% --- Executes on button press in scanforward.
function scanforward_Callback(hObject, eventdata, handles)
% hObject    handle to scanforward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

xlimit=get(gca,'xlim');
handles.xlimit=[xlimit(2) 2*xlimit(2)-xlimit(1)];
xlim(handles.xlimit);
guidata(hObject, handles);


% --- Executes on button press in regret.
function regret_Callback(hObject, eventdata, handles)
% hObject    handle to regret (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% 
% subplot(handles.spectrogram);
% 
% handles.motif_no=handles.motif_no-1;
% guidata(hObject, handles);



% --- Executes on button press in cont.
function cont_Callback(hObject, eventdata, handles)
% hObject    handle to cont (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[handles.cellName,handles.cellPath]=uigetfile('*.mat');
load([handles.cellPath, handles.cellName]);
handles.cell=presentCell;
handles.save=1;
handles.counter=1;
handles.no_syll=size(handles.cell.syll,2)/2;
if (handles.no_syll>1)
    handles.motif_no=length(handles.cell.rec);
    handles.song_no=length(handles.cell.rec);
else
    handles.song_no=length(handles.cell.rec);
end
guidata(hObject, handles);



% --- Executes on button press in delete_pause.
function delete_pause_Callback(hObject, eventdata, handles)
% hObject    handle to delete_pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.counter>1)
    counter=abs(rem(handles.counter-1,2)-2);
    index=find(handles.syllEndTimes < handles.markertime(counter),1,'last');
    handles.syllEndTimes=[handles.syllEndTimes(1:index-1),handles.syllEndTimes(index+1:end)];
    handles.syllStartTimes=[handles.syllStartTimes(1:index),handles.syllStartTimes(index+2:end)];
    delete(handles.linestart(index+1));delete(handles.lineend(index));
    handles.linestart=[handles.linestart(1:index),handles.linestart(index+2:end)];
    handles.lineend=[handles.lineend(1:index-1),handles.lineend(index+1:end)];
end
guidata(hObject, handles);



% --- Executes on button press in AddSong.
 function AddSong_Callback(hObject, eventdata, handles)
% hObject    handle to AddSong (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (length(handles.markertime)>1 & handles.save==1 & handles.spikes==1)
    button='Go Ahead';
    handles.song_no=handles.song_no+1;
    ncounter=abs(rem(handles.counter-1,2)-2)
    delete(handles.marker(1)); 
    delete(handles.marker(2));
    startSong=min(handles.markertime(1),handles.markertime(2))
    endSong=max(handles.markertime(1),handles.markertime(2))
    handles.cell.ID=handles.birdname;
    handles.cell.rec(handles.song_no)=handles.rec;
    handles.cell.motif(handles.song_no,1)=startSong;
    handles.cell.motif(handles.song_no,2)=endSong;
    no_syll=1;
    
    if (handles.motif_no>1 & no_syll~=handles.no_syll)
        button = questdlg('The number of syllables is different from previous motifs','Warning!','New','Go Ahead','New');
    end
    if (handles.motif_no>1 & strcmp(button, 'New'))
            handles.song_no=handles.song_no-1;
    else
    
       handles.no_syll=no_syll;
       handles.cell.syll(handles.song_no,1)=startSong;
       handles.cell.syll(handles.song_no,2)=endSong;
        
        if (handles.spikes==1)
            spike_index=find(handles.sorted_spikes>startSong & handles.sorted_spikes<endSong);
            handles.cell.spikes(handles.song_no,1:(spike_index(end)-spike_index(1))+1)=handles.sorted_spikes(spike_index(1):spike_index(end));
        end
        presentCell=handles.cell;
        save ([handles.cellPath,handles.cellName], 'presentCell');
        %add a screen pointer to the motif that has already been selected
        subplot(handles.spectrogram);
        handles.song_line=line([startSong endSong],[6000 6000],'color','r');
        handles.song_text=text((startSong+endSong)/2,8000,['Song Number ' num2str(handles.song_no)],'HorizontalAlignment','center','Color',[1 1 0]);

    end
    handles.counter=1
    guidata(hObject, handles);
else
    warndlg('Load spikes (or cell) to add Motif')
end


