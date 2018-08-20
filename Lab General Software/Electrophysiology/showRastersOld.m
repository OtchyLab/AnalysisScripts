function varargout = showRastersOld(varargin)
% SHOWRASTERSOLD M-file for showRastersOld.fig
%      SHOWRASTERSOLD, by itself, creates a new SHOWRASTERSOLD or raises the existing
%      singleton*.
%
%      H = SHOWRASTERSOLD returns the handle to a new SHOWRASTERSOLD or the handle to
%      the existing singleton*.
%
%      SHOWRASTERSOLD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOWRASTERSOLD.M with the given input arguments.
%
%      SHOWRASTERSOLD('Property','Value',...) creates a new SHOWRASTERSOLD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before showRastersOld_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to showRastersOld_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help showRastersOld

% Last Modified by GUIDE v2.5 20-Mar-2007 14:31:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @showRastersOld_OpeningFcn, ...
                   'gui_OutputFcn',  @showRastersOld_OutputFcn, ...
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


% --- Executes just before showRastersOld is made visible.
function showRastersOld_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to showRastersOld (see VARARGIN)

% Choose default command line output for showRastersOld
handles.output = hObject;
handles.display=0;
handles.corr_index=0;
handles.corr=[];
handles.ifiring=0;
handles.raw=0;
handles.first_time=0;
handles.local=0;
handles.correlation=0;
handles.in_index=0;
handles.out_index=0;
set(handles.timewarping,'Value',1);
linkaxes([handles.raster handles.spike_stats],'x');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes showRastersOld wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = showRastersOld_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadCell.
function loadCell_Callback(hObject, eventdata, handles)
% hObject    handle to loadCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.display=0;
handles.cell=[0];
handles.bin=[];
handles.xscale=[];
handles.yscale=[];
handles.inst_firing=[];
handles.spikes=[];
[handles.filename,handles.dirpath]=uigetfile('*.mat');
load([handles.dirpath handles.filename]);
cd (handles.dirpath);

handles.cell=presentCell;
guidata(hObject, handles);

% --- Executes on button press in local_linear.
function local_linear_Callback(hObject, eventdata, handles)
% hObject    handle to local_linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (get(hObject,'Value') == get(hObject,'Max'))
    handles.display=2;
else
    handles.display=1;
end
handles.no_syll=size(handles.cell.syll,2)/2;
origin=ceil(handles.no_syll/2)*2-1;
syll_times=handles.cell.syll;
for i=1:2*handles.no_syll
    syll_times(:,i)=handles.cell.syll(:,i)-handles.cell.syll(:,origin);
    handles.motif_template(i)=mean(syll_times(:,i));
end
subplot(handles.raster);
cla;
if (handles.first_time==0)
    handles.first_time=1;
    xlim([(handles.motif_template(1)-0.1) (handles.motif_template(end)+0.1)]);
end 
if (handles.corr_index==1 | handles.ifiring==1)
    xlim([(handles.motif_template(1)-0.1) (handles.motif_template(end)+0.1)]);
    handles.corr_index=0;
    handles.ifiring=0;
end
handles.spikes=cell(1,handles.to_index-handles.from_index+1);
ylim([handles.from_index-1, handles.to_index]);
set(gca,'YDir','normal');
for i=handles.from_index:handles.to_index  %remove the zeros, then subtract the reference point  
    raw_raster(i,1:size(handles.cell.spikes,2))=handles.cell.spikes(i,1:size(handles.cell.spikes,2));
    handles.spikes{i}=nonzeros(raw_raster(i,:))';
    handles.spikes{i}=handles.spikes{i}-handles.cell.syll(i,origin);
    handles.spike_train=handles.spikes{i};
    %startwarping at this marker, should be 1 under normal circumstances
    if (get(handles.timewarping,'Value')==1)
        marker_times=syll_times(i,1:2:end);
        motif_template=handles.motif_template(1:2:end);
        handles.spike_train = warpSpikeTimes(handles.spike_train, marker_times, motif_template);
        handles.spikes{i}=handles.spike_train;
        %also do spike warping
%         for j=1:length(handles.spike_train)
%             line([handles.spike_train(j) handles.spike_train(j)], [i-0.9 i-0.1]);
%         end
    else
        for j=1:length(handles.spike_train)
        line([handles.spike_train(j) handles.spike_train(j)], [i-0.9 i-0.1]);
        end
    end
    
end    
handles.xlimit=get(gca,'xlim');
handles.ylimit=get(gca,'ylim');

%plot the psth
handles.bin=[];
handles.bin_size=10; %bin size for psth in ms
handles.motif_length=(handles.motif_template(end)-handles.motif_template(1)+0.05)*1000; %length in ms
for i=1:floor(handles.motif_length/handles.bin_size)
    handles.bin(i)=handles.motif_template(1)-0.05+(i*handles.bin_size)/1000;
end
handles.bin_size=handles.bin(2)-handles.bin(1);
if (get(handles.timewarping,'Value')==1)
    spikes=cat(2,handles.spikes{handles.from_index:handles.to_index}); 
    psth = histc(spikes,handles.bin)/((handles.to_index-handles.from_index+1)*handles.bin_size);
    startTime=handles.motif_template(1)-0.04;
    endTime=handles.motif_template(end);
    sampleRate=40000;
    GaussWidth=0.008; %in seconds
    motifShifts = crossCorrelate(handles.from_index,handles.to_index, handles.spikes,sampleRate, 0, 0, GaussWidth, startTime,endTime,floor(1/handles.bin_size),psth, 1);
    subplot(handles.raster);
    cla;
    xlim([(handles.motif_template(1)-0.1) (handles.motif_template(end)+0.1)]);
    ylim([handles.from_index-1,handles.to_index]);
    set(gca,'YDir','normal');

    for i=handles.from_index:handles.to_index
        handles.spikes{i}=handles.spikes{i}+ motifShifts(i)/1000;
    end

    for i=handles.from_index:handles.to_index
        from=i-0.9;
        to=i-0.1;
        spikeIndex=find(handles.spikes{i}>=(handles.motif_template(1)-0.1) & handles.spikes{i}<=(handles.motif_template(end)+0.1));
        spikes=handles.spikes{i}(spikeIndex);

        for k=1:length(spikes)
            line([spikes(k) spikes(k)], [from to]);
        end
    end
end
  



if (handles.in_index~=0 & handles.out_index~=0)
    subplot(handles.raster);
    handles.in_line=line(handles.xlimit,[handles.in_index-1 handles.in_index-1],'color','r');
    handles.out_line=line(handles.xlimit,[handles.out_index handles.out_index],'color','r');
    for i=handles.in_index:handles.out_index
        spike_train=handles.spikes{i};
        psth = histc(spike_train,handles.bin);
        if (i==handles.in_index)
            handles.psth_drug=psth;
        else
            handles.psth_drug=handles.psth_drug+psth;
        end
    end
    handles.psth_drug=handles.psth_drug/((handles.out_index-handles.in_index+1)*(handles.bin_size));
    subplot(handles.spike_stats);
    plot(handles.bin,handles.psth_drug,'color','r');

    for i=handles.from_index:handles.in_index-1
        spike_train=handles.spikes{i};
        psth = histc(spike_train,handles.bin);
        if (i==handles.from_index)
            handles.psth_nodrug=psth;
        else
            handles.psth_nodrug=handles.psth_nodrug+psth;
        end
    end
    if (handles.out_index<handles.to_index)
        for i=handles.out_index+1:handles.to_index
            spike_train=handles.spikes{i};
            psth = histc(spike_train,handles.bin);
            handles.psth_nodrug=handles.psth_nodrug+psth;
        end
    end
    handles.psth_nodrug=handles.psth_nodrug/(((handles.in_index-handles.from_index)+(handles.to_index-handles.out_index+1))*(handles.bin_size));
    subplot(handles.spike_stats);
    hold on;
    plot(handles.bin,handles.psth_nodrug);
    hold off;
else 
    for i=handles.from_index:handles.to_index
        spike_train=handles.spikes{i};
        psth = histc(spike_train,handles.bin);
        if (i==handles.from_index)
            handles.psth_nodrug=psth;
        else
            handles.psth_nodrug=handles.psth_nodrug+psth;
        end
    end
    handles.psth_nodrug=handles.psth_nodrug/((handles.to_index-handles.from_index+1)*(handles.bin_size));
    subplot(handles.spike_stats);
    plot(handles.bin,handles.psth_nodrug);
end
sumFR=sum(handles.psth_nodrug(1:end));
         BigN=length(handles.psth_nodrug(1:end));
        
         
         for i=1:BigN
            p(i)=handles.psth_nodrug(i)/sumFR;
            if p(i)>0
                sparseness(i)=(p(i)*log(p(i)));
            else
                sparseness(i)=0;
            end
         end
         sparseness=1+(sum(sparseness)/log(BigN))

xlim([(handles.motif_template(1)-0.1) (handles.motif_template(end)+0.1)]);

guidata(hObject, handles);



function from_Callback(hObject, eventdata, handles)
% hObject    handle to from (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of from as text
%        str2double(get(hObject,'String')) returns contents of from as a double
handles.from_rec=str2double(get(hObject,'String'));
set(handles.from,'String',num2str(handles.from_rec));
handles.from_index=find(handles.cell.rec>=handles.from_rec,1,'first');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function from_CreateFcn(hObject, eventdata, handles)
% hObject    handle to from (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function to_Callback(hObject, eventdata, handles)
% hObject    handle to to (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of to as text
%        str2double(get(hObject,'String')) returns contents of to as a double

handles.to_rec=str2double(get(hObject,'String'));
set(handles.to,'String',num2str(handles.to_rec));
handles.to_index=find(handles.cell.rec<=handles.to_rec,1,'last');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function to_CreateFcn(hObject, eventdata, handles)
% hObject    handle to to (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






function Drug_in_Callback(hObject, eventdata, handles)
% hObject    handle to Drug_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Drug_in as text
%        str2double(get(hObject,'String')) returns contents of Drug_in as a double

handles.in=str2double(get(hObject,'String'));
set(handles.Drug_in,'String',num2str(handles.in));
subplot(handles.raster);
handles.xlimit=get(gca,'xlim');
handles.in_index=find(handles.cell.rec>=handles.in,1,'first');
if (handles.corr_index==1 & handles.in_index>=handles.from_index & handles.in_index<=handles.to_index)
    handles.in_line(1)=line([(handles.in_index-0.5) (handles.in_index-0.5)],[0 handles.to_index+1]);
    handles.in_line(2)=line([0 handles.to_index+1], [(handles.in_index-0.5) (handles.in_index-0.5)]);
else   
    handles.in_line=line(handles.xlimit,[handles.in_index-1 handles.in_index-1],'color','r');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Drug_in_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Drug_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Wash_out_Callback(hObject, eventdata, handles)
% hObject    handle to Wash_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Wash_out as text
%        str2double(get(hObject,'String')) returns contents of Wash_out as a double
handles.out=str2double(get(hObject,'String'));
set(handles.Wash_out,'String',num2str(handles.out));
subplot(handles.raster);
handles.xlimit=get(gca,'xlim');
handles.out_index=find(handles.cell.rec<=handles.out,1,'last');
if (handles.corr_index==1 & handles.out_index>=handles.from_index & handles.out_index<=handles.to_index)
    handles.out_line=line([(handles.out_index-0.5) (handles.out_index-0.5)],[0 handles.to_index+1]);
    handles.out_line=line([0 handles.to_index+1], [(handles.out_index-0.5) (handles.out_index-0.5)]);
else   
    handles.out_line=line(handles.xlimit,[handles.out_index handles.out_index],'color','r');
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Wash_out_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Wash_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in alignment.
function alignment_Callback(hObject, eventdata, handles)
% hObject    handle to alignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
subplot(handles.raster);
handles.ylimit=get(gca,'ylim');
switch (handles.display)
    case 1
        handles.a_line=line([0 0],handles.ylimit, 'color',[0.2 0 0],'LineWidth',0.1);
    case 2
        for i=1:2*handles.no_syll
            handles.a_line(i)=line([handles.motif_template(i) handles.motif_template(i)], handles.ylimit,'color',[1 0.7 0.7],'LineWidth',0.1);
        end
        
end
        
guidata(hObject, handles);


% --- Executes on button press in Correlation.
function Correlation_Callback(hObject, eventdata, handles)
% hObject    handle to Correlation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.corr_index=1;
startTime=handles.motif_template(1)-0.03;
endTime=handles.motif_template(end)+0.02;
sampleRate=40000;
GaussWidth=0.008; %in seconds
handles.corr = pairWiseCorrelations(handles.from_index,handles.to_index,handles.spikes, sampleRate, 0, 0, GaussWidth, startTime,endTime);
%handles.corr=handles.corr(handles.from_index:handles.to_index,handles.from_index:handles.to_index);
subplot(handles.raster);
% handles.xlimit=get(gca,'xlim');
% handles.ylimit=get(gca,'ylim');
imagesc([handles.from_index handles.to_index],[handles.from_index handles.to_index],handles.corr);
set(gca,'YDir','normal');
%colormap('gray');
%write the average correlation.
sum_corr=sum(sum(handles.corr))-sum(diag(handles.corr));
[m,n] = size(handles.corr);
correlation=sum_corr/((m*n)-m)
%make a plot of the sliding correlation
window=2;
for i=window+1:m-window-1
    sliding_corr=handles.corr(i-window:i+window,i-window:i+window);% if window is 2 sliding window is 5....
    sum_sliding_corr=sum(sum(sliding_corr))-sum(diag(sliding_corr));
    average_corr(i)=sum_sliding_corr/((2*window)*(2*window+1));
end
average_corr(1:window)=average_corr(window+1);
average_corr(m-window:m)=average_corr(m-window-1);
subplot(handles.spike_stats);
x=handles.from_index:1:handles.to_index;
plot(x,average_corr,'color','r');
guidata(hObject, handles);

% --- Executes on button press in save_corr.
function save_corr_Callback(hObject, eventdata, handles)
% hObject    handle to save_corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.corrName,handles.corrPath]=uiputfile('corr.mat','Save correlation as:');
presentCorr=handles.corr;
index=[handles.in_index handles.out_index];
save ([handles.corrPath,handles.corrName], 'presentCorr','index');
guidata(hObject, handles);



% --- Executes on button press in inst_firing.
function inst_firing_Callback(hObject, eventdata, handles)
% hObject    handle to inst_firing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.display~=0)
    handles.ifiring=1;
    handles.motif_length=(handles.motif_template(end)-handles.motif_template(1)+0.2)*1000; %length in ms
    handles.no_motifs=(handles.to_index-handles.from_index)+1;
    handles.inst_firing = zeros(handles.no_motifs,handles.motif_length);
    for i=1:handles.no_motifs
        spike_train=handles.spikes{handles.from_index+(i-1)};
        for j=1:handles.motif_length
            if (handles.motif_template(1)-0.1+(j/1000)>=max(spike_train) | handles.motif_template(1)-0.1+(j/1000)<=min(spike_train))
               handles.inst_firing(i,j)=0;
            else
                spike_index=find(spike_train>(handles.motif_template(1)-0.1+(j/1000)),1,'first');
                handles.inst_firing(i,j)=1/(spike_train(spike_index)-spike_train(spike_index-1));
        
            end
         handles.xscale(j)=handles.motif_template(1)-0.1+(j/1000);
        end
        handles.yscale(i)=handles.from_index+i-0.5;
    end
    subplot(handles.raster);
    maxcolor=min([400 max(max(handles.inst_firing))]);
    imagesc(handles.xscale,handles.yscale,handles.inst_firing,[0 maxcolor]);
    set(gca,'YDir','normal');
    guidata(hObject, handles);
else
    warndlg('First plot the spike trains!')
end


% --- Executes on button press in victor.
function victor_Callback(hObject, eventdata, handles)
% hObject    handle to victor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% for i=handles.from_index:handles.to_index
handles.ifiring=1;
handles.no_motif=handles.to_index-handles.from_index+1;
victor=zeros(handles.no_motif,handles.no_motif);
for i=1:handles.no_motif
    spikes_1=handles.spikes{i+handles.from_index-1};
    song_spikes_1=spikes_1(find(handles.motif_template(1)-0.05 < spikes_1 & spikes_1 < handles.motif_template(end)+0.05));
    for j=1:handles.no_motif
        spikes_2=handles.spikes{j+handles.from_index-1};
        song_spikes_2=spikes_2(find(handles.motif_template(1)-0.05 < spikes_2 & spikes_2 < handles.motif_template(end)+0.05));
        victor(i,j)=spkd(song_spikes_1,song_spikes_2,20);
    end
end
subplot(handles.raster);
imagesc([handles.from_index handles.to_index],[handles.from_index handles.to_index],victor);
set(gca,'YDir','normal');
window=2;
[m,n] = size(victor);
sum_victor=sum(sum(victor))-sum(diag(victor));
victor_sum=sum_victor/((m*n)-m)
for i=window+1:m-window-1
    sliding_victor=victor(i-window:i+window,i-window:i+window);% if window is 2 sliding window is 5....
    sum_sliding_victor=sum(sum(sliding_victor))-sum(diag(sliding_victor));
    average_victor(i)=sum_sliding_victor/((2*window)*(2*window+1));
end
average_victor(1:window)=average_victor(window+1);
average_victor(m-window:m)=average_victor(m-window-1);
subplot(handles.spike_stats);
x=handles.from_index:1:handles.to_index;
plot(x,average_victor,'color','r');
guidata(hObject, handles);


% --- Executes on button press in meansubtracted_IF.
function meansubtracted_IF_Callback(hObject, eventdata, handles)
% hObject    handle to meansubtracted_IF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.ifiring==1)
    window=4;
    IF_core(1:handles.out_index-handles.in_index+1,:)=handles.inst_firing(handles.in_index-handles.from_index+1:handles.out_index-handles.from_index+1,:);
    IF_core_mean=sum(IF_core,1)/(handles.out_index-handles.in_index+1);
    for i=ceil(window/2)+1:handles.motif_length-ceil(window/2)%sliding average
        IF_core_mean(i)=sum(IF_core_mean(i-floor(window/2):i+floor(window/2)))/(floor(window/2)*2+1);
    end
    for i=1:handles.no_motifs
        IF=handles.inst_firing(i,:);
        meansubtracted_IF(i,:)=IF-IF_core_mean;
    end
    subplot(handles.raster);
    maxcolor=min([400 max(max(meansubtracted_IF))]);
    mincolor=max([-400 min(min(meansubtracted_IF))]);
    imagesc(handles.xscale,handles.yscale,meansubtracted_IF,[mincolor maxcolor]);
    set(gca,'YDir','normal');
    guidata(hObject, handles);
else
    
    warndlg('First plot the instanteneous firing rate!')
end



% --- Executes on button press in extra_spikes.
function extra_spikes_Callback(hObject, eventdata, handles)
% hObject    handle to extra_spikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

subplot(handles.raster);
cla;
for i=handles.from_index:handles.to_index  %remove the zeros, then subtract the reference point  
   
    handles.spike_train=handles.spikes{i};
    
    for j=1:length(handles.no_motif)
        index=find(handles.bin<=(handles.spike_train(j)),1,'last');
        if (handles.psth_drug(index)==0)
            line([handles.spike_train(j) handles.spike_train(j)], [i-0.9 i-0.1],'color','r');
        end
    end
end   
handles.xlimit=get(gca,'xlim');
handles.ylimit=get(gca,'ylim');
guidata(hObject, handles);

% --- Executes on button press in timewarping.
function timewarping_Callback(hObject, eventdata, handles)
% hObject    handle to timewarping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of timewarping


guidata(hObject, handles);


function saveSound_Callback(hObject, eventdata, handles)
% hObject    handle to saveSound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fs=40000;
motif='';
[file,path] = uigetfile('*.mat', 'Choose an experiment file to load:');
if isequal(file,0) || isequal(path,0)
else
    experName= [path,file];
    load(experName);
    handles.exper=exper;
end  

%for numSyll=1:1%length(handles.sequence)
    old_dir=cd;
    [file,path]=uiputfile('',['Choose folder to save motif ',motif],'motif');
    cd (path);
    for i=1:handles.to_index-handles.from_index+1
        
        sound = loadAudio(exper, handles.cell.rec(handles.from_index+i-1));
        cd(handles.dirpath);
        startTime=handles.cell.motif(handles.from_index+i-1,1)-0.01;
        endTime=handles.cell.motif(handles.from_index+i-1,2)+0.01;
        sylltosave=sound(round(startTime*handles.fs):round(endTime*handles.fs));
        y = resample(sylltosave,44100,handles.fs);
        cd (path);
        if (i>9)
            add='9';
        else
            add='';
        end
        wavwrite(y,44100,[num2str(handles.cell.rec(handles.from_index+i-1)),'_motif',motif,'_',add,num2str(i)]);
        cd(handles.dirpath);

    end