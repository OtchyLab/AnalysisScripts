function varargout = showRasters(varargin)
% SHOWRASTERS M-file for showRasters.fig
%      SHOWRASTERS, by itself, creates a new SHOWRASTERS or raises the existing
%      singleton*.
%
%      H = SHOWRASTERS returns the handle to a new SHOWRASTERS or the handle to
%      the existing singleton*.
%
%      SHOWRASTERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOWRASTERS.M with the given input arguments.
%
%      SHOWRASTERS('Property','Value',...) creates a new SHOWRASTERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before showRasters_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to showRasters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help showRasters

% Last Modified by GUIDE v2.5 04-Dec-2006 15:03:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @showRasters_OpeningFcn, ...
                   'gui_OutputFcn',  @showRasters_OutputFcn, ...
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


% --- Executes just before showRasters is made visible.
function showRasters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to showRasters (see VARARGIN)

% Choose default command line output for showRasters

set(handles.timewarping,'Value',1);
linkaxes([handles.raster handles.spike_stats handles.specgram],'x');
set(handles.specgram,'XTick',[]);
set(handles.specgram,'YTick',[]);
set(handles.popupSequence,'Value',1);
% Update handles structure
handles.corr_index=0;
handles.ifiring=0;
handles.first_time=0;
handles.in_index=0;
handles.out_index=0;
guidata(hObject, handles);

% UIWAIT makes showRasters wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = showRasters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output; %not sure what this does, but i commented
%it out b/c it gave me an error(Bence)


% --- Executes on button press in loadCell.
function loadCell_Callback(hObject, eventdata, handles)
% hObject    handle to loadCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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
if (~isfield(handles,'chosenStartSeq'))
    warndlg('chose the sequence to display first');
    return;
else
    if (get(handles.timewarping,'Value') == get(hObject,'Max'))
        handles.display=2;
    else
        handles.display=1;
    end
    [handles.spikes handles.motif handles.syll]=fillUpRaster(handles);%gives back a spike array, and syllable times
    origin=floor(handles.numberSyll/2)+1;
    handles.noMotifs=size(handles.motif,2);
    for i=1:handles.noMotifs
        syll_times(:,i)=handles.syll(:,i)-handles.syll(origin,i);
    end
    handles.motif_template=[];
    for i=1:handles.numberSyll*2
        handles.motif_template(i)=mean(syll_times(i,:));
    end
    subplot(handles.raster);
    cla;
    xlim([(handles.motif_template(1)-0.1) (handles.motif_template(end)+0.1)]);
    ylim([0,handles.noMotifs]);
    set(gca,'YDir','normal');
    for i=1:handles.noMotifs    
        handles.spikes{i}=handles.spikes{i}-handles.syll(origin,i);
        handles.spike_train=handles.spikes{i};
        %startwarping at this marker, should be 1 under normal circumstances
        if (get(handles.timewarping,'Value')==1)
            marker_times=syll_times(1:2:end,i);
            motif_template=handles.motif_template(1:2:end);
            handles.spike_train = warpSpikeTimes(handles.spike_train, marker_times, motif_template);
            handles.spikes{i}=handles.spike_train; 
            for j=1:length(handles.spike_train)
                line([handles.spike_train(j) handles.spike_train(j)], [i-0.9 i-0.1]);
            end
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
     bin_size=3; %bin size for psth in ms
     motif_length=(handles.motif_template(end)-handles.motif_template(1)+0.8)*1000; %length in ms
     for i=1:floor(motif_length/bin_size)
         handles.bin(i)=handles.motif_template(1)-0.4+(i*bin_size)/1000;
     end
     bin_size=handles.bin(2)-handles.bin(1);
%     if (handles.in_index~=0 & handles.out_index~=0)
%         subplot(handles.raster);
%         handles.in_line=line(handles.xlimit,[handles.in_index-1 handles.in_index-1],'color','r');
%         handles.out_line=line(handles.xlimit,[handles.out_index handles.out_index],'color','r');
%         for i=handles.in_index:handles.out_index
%             spike_train=handles.spikes{i};
%             psth = histc(spike_train,handles.bin);
%             if (i==handles.in_index)
%                 handles.psth_drug=psth;
%             else
%                 handles.psth_drug=handles.psth_drug+psth;
%             end
%         end
%         handles.psth_drug=handles.psth_drug/((handles.out_index-handles.in_index+1)*(handles.bin_size));
%         subplot(handles.spike_stats);
%         plot(handles.bin,handles.psth_drug,'color','r');
% 
%         for i=handles.from_index_audio_audio:handles.in_index-1
%             spike_train=handles.spikes{i};
%             psth = hist(spike_train,handles.bin);
%             if (i==handles.from_index_audio_audio)
%                 handles.psth_nodrug=psth;
%             else
%                 handles.psth_nodrug=handles.psth_nodrug+psth;
%             end
%         end
%         if (handles.out_index<handles.to_index_audio)
%             for i=handles.out_index+1:handles.to_index_audio
%                 spike_train=handles.spikes{i};
%                 psth = hist(spike_train,handles.bin);
%                 handles.psth_nodrug=handles.psth_nodrug+psth;
%             end
%         end
%         handles.psth_nodrug=handles.psth_nodrug/(((handles.in_index-handles.from_index_audio_audio)+(handles.to_index_audio-handles.out_index+1))*(handles.bin_size));
%         subplot(handles.spike_stats);
%         hold on;
%         plot(handles.bin,handles.psth_nodrug);
%         hold off;
%     else 
         for i=1:handles.noMotifs
             spike_train=handles.spikes{i};
             psth = hist(spike_train,handles.bin);
             if (i==1)
                 handles.psth_nodrug=psth;
             else
                 handles.psth_nodrug=handles.psth_nodrug+psth;
             end
         end
         handles.psth_nodrug=handles.psth_nodrug/(handles.noMotifs*bin_size);
         subplot(handles.spike_stats);
         plot(handles.bin,handles.psth_nodrug);
%     end

     xlim([(handles.motif_template(1)-0.1) (handles.motif_template(end)+0.1)]);
end
guidata(hObject, handles);

function [spikes motif syll]=fillUpRaster(handles)

for i=1:size(handles.chosenStartSeq,1)
    syll(1:2:handles.numberSyll*2-1,i)=handles.audioAnnotation{1,handles.chosenStartSeq(i,1)}.segFileStartTimes(handles.chosenStartSeq(i,2):handles.chosenStartSeq(i,2)+handles.numberSyll-1);
    syll(2:2:handles.numberSyll*2,i)=handles.audioAnnotation{1,handles.chosenStartSeq(i,1)}.segFileEndTimes(handles.chosenStartSeq(i,2):handles.chosenStartSeq(i,2)+handles.numberSyll-1);
    motif(1,i)= handles.audioAnnotation{1,handles.chosenStartSeq(i,1)}.segFileStartTimes(handles.chosenStartSeq(i,2));
    motif(2,i)= handles.audioAnnotation{1,handles.chosenStartSeq(i,1)}.segFileEndTimes(handles.chosenStartSeq(i,2)+handles.numberSyll-1);
    
    presentSpikeIndex=find(handles.spikeIndex==handles.audioAnnotation{1,handles.chosenStartSeq(i,1)}.filenum);
    if isempty(presentSpikeIndex)
        warndlg('no spikes for file %d',handles.audioAnnotation{1,handles.chosenStartSeq(i,1)}.filenum);
        spikes{i}=[];
    else
        spikes{i}=handles.cell(1,presentSpikeIndex).spikes;
    end
    
end
    

function from_Callback(hObject, eventdata, handles)
% hObject    handle to from (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of from as text
%        str2double(get(hObject,'String')) returns contents of from as a double
if (isfield(handles, 'filenums'))
    handles.from_rec=str2double(get(hObject,'String'));
    set(handles.from,'String',num2str(handles.from_rec));
    handles.from_index_audio=find(handles.filenums>=handles.from_rec,1,'first');
    guidata(hObject, handles);
else
    warndlg('Hey! Load a valid annotation file first')
end

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

if (isfield(handles, 'filenums'))
    handles.to_rec=str2double(get(hObject,'String'));
    set(handles.to,'String',num2str(handles.to_rec));
    handles.to_index_audio=find(handles.filenums>=handles.to_rec,1,'first');
    guidata(hObject, handles);
else
    warndlg('Hey! Load a valid annotation file first')
end


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
        for i=1:2*handles.numberSyll
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
handles.corr = pairWiseCorrelations(1,handles.noMotifs, handles.spikes,sampleRate, 0, 0, GaussWidth, startTime,endTime);
%handles.corr=handles.corr(handles.from_index_audio:handles.to_index,handles.from_index_audio:handles.to_index);
subplot(handles.raster);
% handles.xlimit=get(gca,'xlim');
% handles.ylimit=get(gca,'ylim');
imagesc([1 handles.noMotifs],[1 handles.noMotifs],handles.corr);
set(gca,'YDir','normal');
%colormap('gray');
%write the average correlation.
sum_corr=sum(sum(handles.corr))-sum(diag(handles.corr));
[m,n] = size(handles.corr);
correlation=sum_corr/((m*n)-m)
%make a plot of the sliding correlation
window=2;
if (m>window*2+1)
    
    for i=window+1:m-window-1
        sliding_corr=handles.corr(i-window:i+window,i-window:i+window);% if window is 2 sliding window is 5....
        sum_sliding_corr=sum(sum(sliding_corr))-sum(diag(sliding_corr));
        average_corr(i)=sum_sliding_corr/((2*window)*(2*window+1));
    end
else
    warndlg('too few motifs to calculate sliding correlation');
    return;
end    
average_corr(1:window)=average_corr(window+1);
average_corr(m-window:m)=average_corr(m-window-1);
subplot(handles.spike_stats);
x=1:handles.noMotifs;
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
    motif_length=(handles.motif_template(end)-handles.motif_template(1)+0.2)*1000; %length in ms
    handles.inst_firing = zeros(handles.noMotifs,motif_length);
    for i=1:handles.noMotifs
        spike_train=handles.spikes{i};
        for j=1:motif_length
            if (handles.motif_template(1)-0.1+(j/1000)>=max(spike_train) | handles.motif_template(1)-0.1+(j/1000)<=min(spike_train))
               handles.inst_firing(i,j)=0;
            else
                spike_index=find(spike_train>(handles.motif_template(1)-0.1+(j/1000)),1,'first');
                handles.inst_firing(i,j)=1/(spike_train(spike_index)-spike_train(spike_index-1));
        
            end
        handles.xscale(j)=handles.motif_template(1)-0.1+(j/1000);
        end
        handles.yscale(i)=i-0.5;
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

% for i=handles.from_index_audio:handles.to_index
handles.ifiring=1;
victor=zeros(handles.noMotifs,handles.noMotifs);
for i=1:handles.noMotifs
    spikes_1=handles.spikes{i};
    song_spikes_1=spikes_1(find(handles.motif_template(1)-0.05 < spikes_1 & spikes_1 < handles.motif_template(end)+0.05));
    for j=1:handles.noMotifs
        spikes_2=handles.spikes{j};
        song_spikes_2=spikes_2(find(handles.motif_template(1)-0.05 < spikes_2 & spikes_2 < handles.motif_template(end)+0.05));
        victor(i,j)=spkd(song_spikes_1,song_spikes_2,20);
    end
end
subplot(handles.raster);
imagesc([1 handles.noMotifs],[1 handles.noMotifs],victor);
set(gca,'YDir','normal');
window=2;
[m,n] = size(victor);
for i=window+1:m-window-1
    sliding_victor=victor(i-window:i+window,i-window:i+window);% if window is 2 sliding window is 5....
    sum_sliding_victor=sum(sum(sliding_victor))-sum(diag(sliding_victor));
    average_victor(i)=sum_sliding_victor/((2*window)*(2*window+1));
end
average_victor(1:window)=average_victor(window+1);
average_victor(m-window:m)=average_victor(m-window-1);
subplot(handles.spike_stats);
x=1:handles.noMotifs;
plot(x,average_victor,'color','r');
guidata(hObject, handles);


% --- Executes on button press in meansubtracted_IF.
function meansubtracted_IF_Callback(hObject, eventdata, handles)
% hObject    handle to meansubtracted_IF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.ifiring==1)
    window=4;
    motif_length=(handles.motif_template(end)-handles.motif_template(1)+0.2)*1000; %length in ms
    IF_core(1:handles.noMotifs,:)=handles.inst_firing(1:handles.noMotifs,:);
    IF_core_mean=sum(IF_core,1)/(handles.noMotifs);
    for i=ceil(window/2)+1:motif_length-ceil(window/2)%sliding average
        IF_core_mean(i)=sum(IF_core_mean(i-floor(window/2):i+floor(window/2)))/(floor(window/2)*2+1);
    end
    for i=1:handles.noMotifs
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

% subplot(handles.raster);
% cla;
% for i=1:handles.noMotifs  %remove the zeros, then subtract the reference point  
%    
%     handles.spike_train=handles.spikes{i};
%     
%     for j=1:length(handles.spike_train)
%         index=find(handles.bin<=(handles.spike_train(j)),1,'last');
%         if (handles.psth_drug(index)==0)
%             line([handles.spike_train(j) handles.spike_train(j)], [i-0.9 i-0.1],'color','r');
%         end
%     end
% end   
% handles.xlimit=get(gca,'xlim');
% handles.ylimit=get(gca,'ylim');
% guidata(hObject, handles);

% --- Executes on button press in timewarping.
function timewarping_Callback(hObject, eventdata, handles)
% hObject    handle to timewarping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of timewarping
guidata(hObject, handles);


% --- Executes on button press in loadAnnotation.
function loadAnnotation_Callback(hObject, eventdata, handles)
% hObject    handle to loadAnnotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%cd (handles.dirpath);
[file,path] = uigetfile('*.mat', 'Choose an audio annotation .mat file to load:');

if isequal(file,0) || isequal(path,0)
else
    %if they don't hit cancel..
    handles.annotationFileName = [path,filesep,file];
    audioAnnotation = aaLoadHashtable(handles.annotationFileName);
    handles.audioAnnotation = audioAnnotation.elements;
    handles.filenums=getFieldVector(handles.audioAnnotation,'filenum');
end
guidata(hObject, handles);


% --- Executes on selection change in popupNumberSyll.
function popupNumberSyll_Callback(hObject, eventdata, handles)
% hObject    handle to popupNumberSyll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupNumberSyll contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupNumberSyll

handles.numberSyll = get(handles.popupNumberSyll, 'Value');
if (~isfield(handles,'to_index_audio') || ~isfield(handles,'from_index_audio'))
    warndlg('set the records that you want to display');
else
    handles=sortSequences(handles); %sort the sequences of handles.numSyll in terms of how often they appear
    histx=(1:10^handles.numberSyll); %the histogram of the sequences.
    hist_seq=hist(handles.allSeq,histx);
    %choose the 5 most abundant sequences to display in the sequence pop-up
    %window
    [sortSeq, sortIndex]=sort(hist_seq,'descend');
    numSeq=length(find(sortSeq>0));
    if isempty(numSeq)
        warndlg('no sequences');
        return;
    elseif (numSeq<5)
        popupNum=numSeq;
    else
        popupNum=5; %display only the four most frequent sequences
    end   
    popupSeqString(1,:)=[num2str(sortIndex(1)),' (n=', ' ',num2str(sortSeq(1)),')'];
    handles.seq=sortIndex(1);
    if (popupNum>1)
        for i=2:popupNum %display the frequency of the various sequences)
            padLength=length(num2str(sortSeq(1)))-length(num2str(sortSeq(i)));
            if (padLength==1)
                freq(i,:)=['  ',num2str(sortSeq(i))];
            elseif (padLength==2)
                freq(i,:)=['   ',num2str(sortSeq(i))];
            else
                freq(i,:)=[' ',num2str(sortSeq(i))];
            end
            handles.seq(i)=sortIndex(i);
            popupSeqString(i,:)=[num2Str(sortIndex(i)),' (n=', freq(i,:),')'];
        end
    end
        
    %set(handles.popupSequence,'String',[num2Str(sortIndex(1)),' (n=', freq(1,:),')';num2Str(sortIndex(2)),' (n=', freq(2,:),')';num2Str(sortIndex(3)),' (n=', freq(3,:),')';num2Str(sortIndex(4)),' (n=', freq(4,:),')']);
    set(handles.popupSequence,'String',popupSeqString);
    set(handles.popupSequence,'Value',1);
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function popupNumberSyll_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupNumberSyll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles=sortSequences(handles)
%make a cell array with all the song sequences

handles.startSeq=[];
[allSeq handles.startSeq counter]=getSequences(handles,handles.numberSyll);
%go through this cell array and build up a array of sequences
%numberSyll long
for i=1:handles.numberSyll
    multiplier(i)=10^(handles.numberSyll-i);
end
handles.allSeq=[];
for i=1:counter-1
    handles.allSeq(i)=dot(allSeq(i,:),multiplier);
end

function [allSeq startSeq counter]=getSequences(handles,numberSyll)
%gets all the sequences containing numSyll
counter=1;
for i=handles.from_index_audio:handles.to_index_audio
    for j=1:length(handles.audioAnnotation{1,i}.segFileStartTimes)-1
        pauses(j)=handles.audioAnnotation{1,i}.segFileStartTimes(j+1)-handles.audioAnnotation{1,i}.segFileEndTimes(j);
        song_interruptions=find(pauses>0.1); %gives the syllable index after which there is a pause (song interruption, no sound for 100ms)

    end
    if (length(handles.audioAnnotation{1,i}.segFileStartTimes)>=numberSyll)
        for j=1:length(handles.audioAnnotation{1,i}.segFileStartTimes)-numberSyll+1
            current_seq(1:numberSyll)=handles.audioAnnotation{1,i}.segType(j:j+numberSyll-1);
            if (isempty(find(current_seq<1 | current_seq>9)) && isempty(find(song_interruptions>j-1 & song_interruptions<j+numberSyll-1)))
                allSeq(counter,1:numberSyll)=current_seq;
                startSeq(counter,:)=[i,j];%index (index in annotation file,number of starting syllable of sequences)
                counter=counter+1;
            end
        end  
    end
end


% --- Executes on selection change in popupSequence.
function popupSequence_Callback(hObject, eventdata, handles)
% hObject    handle to popupSequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupSequence contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupSequence

chosenSeq = handles.seq(get(handles.popupSequence, 'Value')); %chosen sequence as number
startIndex=find(handles.allSeq==chosenSeq); %finds the index of the chosen sequence as they appear in handles.all_seq
handles.chosenStartSeq=[];
for i=1:length(startIndex)
    handles.chosenStartSeq(i,:)=handles.startSeq(startIndex(i),:); 
end
%handles.ChosenStartSeq is a (n,2) matrix where n is the number of
%instances for the chosen sequence, and (i,1) is the filenumber for the ith
%seq, (i,2) the syllable in that file that indicates the start of the ith
%seqence
%fill up the syllable and motif matrices 
handles.spikeIndex=getFieldVector(handles.cell,'filenum');%gets the filenumbers that go with the index in the cell files
sequence(1:handles.numberSyll)=handles.audioAnnotation{1,handles.chosenStartSeq(1,1)}.segFileType(handles.chosenStartSeq(1,2):handles.chosenStartSeq(1,2)+handles.numberSyll)
%displayTemplate(handles,sequence);
guidata(hObject, handles);



function displayTemplate(handles);


% --- Executes during object creation, after setting all properties.
function popupSequence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupSequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SequenceMatrix.
function SequenceMatrix_Callback(hObject, eventdata, handles)
% hObject    handle to SequenceMatrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%go through the annotations and whenever there is a stretch of contiguous
%sequence assign the transition to the matrix. display the figure.
numberSyll=2;
[allSeq startSeq counter]=getSequences(handles,numberSyll);
numberSyllTypes=max(allSeq);
seqMatrix=zeros(numberSyllTypes,numberSyllTypes);   
for i=1:counter-1
    from=allSeq(i,1);
    to=allSeq(i,2);
    seqMatrix(from,to)=seqMatrix(from,to)+1;
end
figure;imagesc([1 numberSyllTypes],[1 numberSyllTypes],seqMatrix);



% --- Executes on button press in loadTemplate.
function loadTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to loadTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.filename,handles.dirpath]=uigetfile('*.mat');
load([handles.dirpath handles.filename]);
cd (handles.dirpath);
handles.templates=templates;
guidata(hObject, handles);

