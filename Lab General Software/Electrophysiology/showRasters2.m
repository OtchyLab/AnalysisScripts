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

% Last Modified by GUIDE v2.5 09-Dec-2006 15:23:46

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
linkaxes([handles.raster handles.spikeStats handles.specgram],'x');
set(handles.specgram,'XTick',[]);
set(handles.specgram,'YTick',[]);
set(handles.raster,'XTick',[]);
set(handles.raster,'YTick',[]);
set(handles.popupSequence,'Value',1);
% Update handles structure
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
if isequal(handles.filename,0) || isequal(handles.dirpath,0)
else
    if isfield(handles,'cell') %if there already exists a cell, add to it or replace it?
        button = questdlg('Do you want to add this cells to previously loaded cells (simultaneously recorded)?','Add a cell','Add','Replace','Add');
        if (strcmp(button,'Replace'))
            handles.cell={};
            load([handles.dirpath handles.filename]);
            cd (handles.dirpath);
            handles.cell{1}=presentCell;
        elseif (strcmp(button,'Add'))
            load([handles.dirpath handles.filename]);
            cd (handles.dirpath);
            handles.cell{length(handles.cell)+1}=presentCell;
        end
    else
        load([handles.dirpath handles.filename]);
        cd (handles.dirpath);
        handles.cell{1}=presentCell;
    end   
end
if isfield(handles,'audioAnnotation')
        
        for i=1:length(handles.cell)    
            from(i)=find(handles.filenums>=handles.cell{i}(1,1).filenum & handles.filenums<=handles.cell{i}(1,end).filenum ,1,'first');
        end
        handles.fromIndexAudio=max(from);
        set(handles.from,'String',num2str(handles.filenums(handles.fromIndexAudio)));
        for i=1:length(handles.cell)    
            to(i)=find(handles.filenums>=handles.cell{i}(1,1).filenum & handles.filenums<=handles.cell{i}(1,end).filenum ,1,'last');
        end
        handles.toIndexAudio=min(to);
        set(handles.to,'String',num2str(handles.filenums(handles.toIndexAudio)));
        audioAnnotation = aaLoadHashtable(handles.annotationFileName);
        handles.audioAnnotation = audioAnnotation.elements;
        handles.filenums=getAnnotationVector(handles.audioAnnotation,'filenum');
end
guidata(hObject, handles);
% --- Executes on button press in displaySpikes.

function displaySpikes_Callback(hObject, eventdata, handles)
% hObject    handle to displaySpikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (~isfield(handles,'chosenStartSeq'))
    warndlg('chose the sequence to display first');
    uiwait;
    return;
else

    if (isfield(handles,'lineOut'))
        handles.lineOut=[];
    end
    if (isfield(handles,'lineIn'))
        handles.lineIn=[];
    end
    color{1}=[1 0 0];
    color{2}=[0 0 0];
    color{3}=[0 0 1];
    color{4}=[0 1 0];
    color{5}=[0 1 1];
    [handles.spikes handles.motif handles.syll]=fillUpRaster(handles);%gives back a spike array, and syllable times
    origin=handles.numberSyll-rem(handles.numberSyll+1,2);
    handles.noMotifs=size(handles.motif,2);
    for i=1:handles.noMotifs
        syllTimes(:,i)=handles.syll(:,i)-handles.syll(origin,i);
        for j=1:length(handles.cell)
            handles.spikes{j,i}=handles.spikes{j,i}-handles.syll(origin,i);
        end
    end
    handles.motifTemplate=[];
    for i=1:handles.numberSyll*2
        handles.motifTemplate(i)=mean(syllTimes(i,:));
    end
    xdata=[handles.motifTemplate(1)-0.1 handles.motifTemplate(end)+0.1];
    displayTemplate(handles,handles.sequence,handles.motifTemplate,xdata)
    subplot(handles.raster);
    cla;
    xlim([(handles.motifTemplate(1)-0.1) (handles.motifTemplate(end)+0.1)]);
    ylim([0,handles.noMotifs]);
    set(gca,'YDir','normal');
    
    for j=1:length(handles.cell)
        for i=1:handles.noMotifs 
            %handles.spikeTrain=handles.spikes{j,i};
            %startwarping at this marker, should be 1 under normal circumstances
            if (get(handles.timewarping,'Value')==1)
                markerTimes=syllTimes(1:2:end,i);
                motifTemplate=handles.motifTemplate(1:2:end);
                handles.spikes{j,i} = warpSpikeTimes(handles.spikes{j,i}, markerTimes, motifTemplate);
                %handles.spikes{j,i}=handles.spikeTrain;
            end
        end
    end
    
    for j=1:length(handles.cell)
        for i=1:handles.noMotifs
            from=(i-1)*length(handles.cell)+(j)-0.9;
            to=(i-1)*length(handles.cell)+j-0.1;
            for k=1:length(handles.spikes{j,i})
                line([handles.spikes{j,i}(k) handles.spikes{j,i}(k)], [from to], 'color', color{j});
            end
        end
    end    
    handles.xlimit=get(gca,'xlim');
    handles.ylimit=get(gca,'ylim');
    if (get(handles.alignment,'Value') == get(handles.alignment,'Max'))
        handles=alignment(handles);
    end
    
    %plot the psth
     handles.bin=[];
     bin_size=3; %bin size for psth in ms
     motif_length=(handles.motifTemplate(end)-handles.motifTemplate(1)+0.8)*1000; %length in ms
     for i=1:floor(motif_length/bin_size)
         handles.bin(i)=handles.motifTemplate(1)-0.4+(i*bin_size)/1000;
     end
     bin_size=handles.bin(2)-handles.bin(1);
     % check wether there has been a drug infusion
     if (size(handles.cell)==1)
         drugStatus=getAnnotationVector(handles.audioAnnotation,'drugindex');
         diffDrugStatus=diff(drugStatus);
         drugIn=find(diffDrugStatus>0)+1; %gives the filenumber at which the switch occurs 
         drug=find(drugStatus==0);
         drug=find(drug>=handles.fromIndexAudio & drug<=handles.toIndexAudio);
         nodrug=find(drugStatus~=0);
         nodrug=find(nodrug>=handles.fromIndexAudio & nodrug<=handles.toIndexAudio);
         drugInIndex=find(drugIn>=handles.fromIndexAudio & drugIn<=handles.toIndexAudio);
         drugIn=drugIn(drugInIndex);
         drugOut=find(diffDrugStatus<0);
         drugOutIndex=find(drugOut>=handles.fromIndexAudio & drugOut<=handles.toIndexAudio);
         drugOut=drugOut(drugOutIndex);
         drug_counter=0;
         nodrug_counter=0;
         if (~isempty(drugIn))
             if (length(drugIn)>length(drugOut))
                 drugOut(end+1)=max(handles.chosenStartSeq(:,1));
             end
             for i=1:length(drugIn)
                handles.lineIn(i)=find(handles.chosenStartSeq(:,1)>=drugIn(i),1,'first');
                handles.lineOut(i)=find(handles.chosenStartSeq(:,1)<=drugOut(i),1,'last');

                if isempty(handles.lineIn)
                    return;
                end
                subplot(handles.raster);
                line(handles.xlimit,[handles.lineIn(i)-1 handles.lineIn(i)-1],'color','r');
                line(handles.xlimit,[handles.lineOut(i) handles.lineOut(i)],'color','k');
             end
         end

         handles.psth_drug=zeros(1,length(handles.bin));
         handles.psth_nodrug=zeros(1,length(handles.bin));
         [c, ai, bi] = intersect((handles.chosenStartSeq(:,1)), nodrug','rows')
         nodrugSpikes=cat(1,handles.spikes{1,ai});
         [c, ai, bi] = intersect((handles.chosenStartSeq(:,1)), drug','rows')
         drugSpikes=cat(1,handles.spikes{1,ai});
         handles.psth_drug = hist(drugSpikes,handles.bin)/(length(drug)*bin_size);
         handles.psth_nodrug = hist(nodrugSpikes,handles.bin)/(length(nodrug)*bin_size);
         subplot(handles.spikeStats);
         cla;
         plot(handles.bin,handles.psth_nodrug,handles.bin,handles.psth_drug,'r');
         xlim([(handles.motifTemplate(1)-0.1) (handles.motifTemplate(end)+0.1)]);
     else
         %if more than 1 cell
         motifs=1:handles.noMotifs;
         subplot(handles.spikeStats);
         cla;
         for j=1:length(handles.cell)
            spikes=cat(1,handles.spikes{j,motifs});
            psth=histc(spikes,handles.bin)/(handles.noMotifs*bin_size);
            plot(handles.bin,psth,'color',color{j});
            xlim([(handles.motifTemplate(1)-0.1) (handles.motifTemplate(end)+0.1)]);
            hold on;
         end
         hold off;
     end
end

guidata(hObject, handles);



function [spikes motif syll]=fillUpRaster(handles)

for i=1:size(handles.chosenStartSeq,1)
    syll(1:2:handles.numberSyll*2-1,i)=handles.audioAnnotation{1,handles.chosenStartSeq(i,1)}.segFileStartTimes(handles.chosenStartSeq(i,2):handles.chosenStartSeq(i,2)+handles.numberSyll-1);
    syll(2:2:handles.numberSyll*2,i)=handles.audioAnnotation{1,handles.chosenStartSeq(i,1)}.segFileEndTimes(handles.chosenStartSeq(i,2):handles.chosenStartSeq(i,2)+handles.numberSyll-1);
    motif(1,i)= handles.audioAnnotation{1,handles.chosenStartSeq(i,1)}.segFileStartTimes(handles.chosenStartSeq(i,2));
    motif(2,i)= handles.audioAnnotation{1,handles.chosenStartSeq(i,1)}.segFileEndTimes(handles.chosenStartSeq(i,2)+handles.numberSyll-1);
    for j=1:length(handles.cell)
        presentSpikeIndex=find(handles.spikeIndex{j}==handles.audioAnnotation{1,handles.chosenStartSeq(i,1)}.filenum);  
        if isempty(presentSpikeIndex)
            warndlg(['no spikes for file ' num2str(handles.audioAnnotation{1,handles.chosenStartSeq(i,1)}.filenum) ' cell number' num2str(j)]);
            spikes{i}=[];
            uiwait;
        else
            spikes{j,i}=handles.cell{j}(1,presentSpikeIndex).spikes;
        end
    end
end
    

function from_Callback(hObject, eventdata, handles)
% hObject    handle to from (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of from as text
%        str2double(get(hObject,'String')) returns contents of from as a double
if (isfield(handles, 'filenums'))
    fromRec=str2double(get(hObject,'String'));
    toRec=str2double(get(handles.to,'String'));
    set(handles.from,'String',num2str(fromRec));
    handles.fromIndexAudio=find(handles.filenums>=fromRec & handles.filenums<=toRec ,1,'first');
    %handles.fromIndexAudio=find(handles.filenums>=fromRec,1,'first');
    guidata(hObject, handles);
else
    warndlg('Hey! Load a valid annotation file first');
    uiwait;
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
    toRec=str2double(get(hObject,'String'));
    fromRec=str2double(get(handles.from,'String'));
    set(handles.to,'String',num2str(toRec));
    handles.toIndexAudio=find(handles.filenums<=toRec & handles.filenums>=fromRec,1,'last');
    guidata(hObject, handles);
else
    warndlg('Hey! Load a valid annotation file first');
    uiwait;
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

guidata(hObject, handles);

function handles=alignment(handles)
subplot(handles.raster);
handles.ylimit=get(gca,'ylim');
if (get(handles.timewarping,'Value') == get(handles.timewarping,'Max'))
     handles.display=2;
else
     handles.display=1;
end
switch (handles.display)
    case 1
        handles.a_line=line([0 0],handles.ylimit, 'color',[0.2 0 0],'LineWidth',0.1);
    case 2
        for i=1:2*handles.numberSyll
            handles.a_line(i)=line([handles.motifTemplate(i) handles.motifTemplate(i)], handles.ylimit,'color',[1 0.7 0.7],'LineWidth',0.1);
        end
        
end


% --- Executes on button press in Correlation.
function Correlation_Callback(hObject, eventdata, handles)
% hObject    handle to Correlation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (isfield(handles,'spikes'))
    if (size(handles.cell)>1)
        warndlg('this function not valid for multiple cells.');
        uiwait;
        return;
    end
    startTime=handles.motifTemplate(1)-0.03;
    endTime=handles.motifTemplate(end)+0.02;
    sampleRate=40000;
    GaussWidth=0.008; %in seconds
    handles.corr = pairWiseCorrelations(1,handles.noMotifs, handles.spikes{1},sampleRate, 0, 0, GaussWidth, startTime,endTime);
    %handles.corr=handles.corr(handles.fromIndexAudio:handles.to_index,handles.fromIndexAudio:handles.to_index);
    subplot(handles.raster);
    % handles.xlimit=get(gca,'xlim');
    % handles.ylimit=get(gca,'ylim');
    imagesc([1 handles.noMotifs],[1 handles.noMotifs],handles.corr);
    set(gca,'YDir','normal');

    %draw drug lines
    if (isfield(handles,'lineIn')&& ~isempty(handles.lineIn)) 
        for i=1:length(handles.lineIn)

            subplot(handles.raster);
            line([(handles.lineIn(i)-0.5) (handles.lineIn(i)-0.5)],[0.5 handles.noMotifs+0.5],'color','r');
            line([0.5 handles.noMotifs+0.5], [(handles.lineIn(i)-0.5) (handles.lineIn(i)-0.5)],'color','r');
            line([(handles.lineOut(i)+0.5) (handles.lineOut(i)+0.5)],[0.5 handles.noMotifs+0.5],'color','k');
            line([0.5 handles.noMotifs+0.5], [(handles.lineOut(i)+0.5) (handles.lineOut(i)+0.5)],'color','k');
    
        end
    end
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
        uiwait;
        return;
    end    
    average_corr(1:window)=average_corr(window+1);
    average_corr(m-window:m)=average_corr(m-window-1);
    subplot(handles.spikeStats);
    cla;
    x=1:handles.noMotifs;
    plot(x,average_corr,'color','r');
    guidata(hObject, handles);
else
    warndlg('plot the spikes first')
    uiwait;
end


% --- Executes on button press in inst_firing.
function inst_firing_Callback(hObject, eventdata, handles)
% hObject    handle to inst_firing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (size(handles.cell)>1)
    warndlg('this function not valid for multiple cells.');
    uiwait;
    return;
end
if (isfield(handles,'motifTemplate'))
    if (size(handles.cell)>1)
        warndlg('this function not valid for multiple cells.');
        uiwait;
        return;
    end
    handles.xscale=[];
    handles.yscale=[];
    motif_length=floor((handles.motifTemplate(end)-handles.motifTemplate(1)+0.2)*1000); %length in ms
    handles.inst_firing = zeros(handles.noMotifs,motif_length);
    for i=1:handles.noMotifs
        spikeTrain=handles.spikes{1,i};
        if ~isempty(spikeTrain)
           
            for j=1:motif_length
                a=(handles.motifTemplate(1)-0.1+(j/1000)>=max(spikeTrain));
                b=(handles.motifTemplate(1)-0.1+(j/1000)<=min(spikeTrain)) ;
                
                if ((a || b))
                   handles.inst_firing(i,j)=0;
                else
                    spike_index=find(spikeTrain>(handles.motifTemplate(1)-0.1+(j/1000)),1,'first');            
                    handles.inst_firing(i,j)=1/(spikeTrain(spike_index)-spikeTrain(spike_index-1));

                end
            handles.xscale(j)=handles.motifTemplate(1)-0.1+(j/1000);
            end
            
        end
        handles.yscale(i)=i-0.5;
    end
    subplot(handles.raster);
    maxcolor=min([400 max(max(handles.inst_firing))]);
    imagesc(handles.xscale,handles.yscale,handles.inst_firing,[0 maxcolor]);
    set(gca,'YDir','normal');
    if (isfield(handles,'lineIn'))
        for i=1:length(handles.lineIn)
            line(handles.xlimit,[handles.lineIn(i)-1 handles.lineIn(i)-1],'color','r');
            line(handles.xlimit,[handles.lineOut(i) handles.lineOut(i)],'color','k');
        end
    end
    if (get(handles.alignment,'Value') == get(handles.alignment,'Max'))
        handles=alignment(handles);
    end
    guidata(hObject, handles);
else
    warndlg('First plot the spike trains!')
    uiwait;
end


% --- Executes on button press in victor.
function victor_Callback(hObject, eventdata, handles)
% hObject    handle to victor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (size(handles.cell)>1)
    warndlg('this function not valid for multiple cells.');
    uiwait;
    return;
end
victor=zeros(handles.noMotifs,handles.noMotifs);
for i=1:handles.noMotifs
    spikes_1=handles.spikes{i};
    song_spikes_1=spikes_1(find(handles.motifTemplate(1)-0.05 < spikes_1 & spikes_1 < handles.motifTemplate(end)+0.05));
    for j=1:handles.noMotifs
        spikes_2=handles.spikes{j};
        song_spikes_2=spikes_2(find(handles.motifTemplate(1)-0.05 < spikes_2 & spikes_2 < handles.motifTemplate(end)+0.05));
        victor(i,j)=spkd(song_spikes_1,song_spikes_2,20);
    end
end
subplot(handles.raster);
imagesc([1 handles.noMotifs],[1 handles.noMotifs],victor);
set(gca,'YDir','normal');
if (isfield(handles,'lineIn')&& ~isempty(handles.lineIn)) 
    for i=1:length(handles.lineIn)
        subplot(handles.raster);
        line([(handles.lineIn(i)-0.5) (handles.lineIn(i)-0.5)],[0.5 handles.noMotifs+0.5],'color','r');
        line([0.5 handles.noMotifs+0.5], [(handles.lineIn(i)-0.5) (handles.lineIn(i)-0.5)],'color','r');
        line([(handles.lineOut(i)+0.5) (handles.lineOut(i)+0.5)],[0.5 handles.noMotifs+0.5],'color','k');
        line([0.5 handles.noMotifs+0.5], [(handles.lineOut(i)+0.5) (handles.lineOut(i)+0.5)],'color','k');
    end
end
window=2;
[m,n] = size(victor);
for i=window+1:m-window-1
    sliding_victor=victor(i-window:i+window,i-window:i+window);% if window is 2 sliding window is 5....
    sum_sliding_victor=sum(sum(sliding_victor))-sum(diag(sliding_victor));
    average_victor(i)=sum_sliding_victor/((2*window)*(2*window+1));
end
average_victor(1:window)=average_victor(window+1);
average_victor(m-window:m)=average_victor(m-window-1);
subplot(handles.spikeStats);
x=1:handles.noMotifs;
plot(x,average_victor,'color','r');
guidata(hObject, handles);


% --- Executes on button press in meansubtracted_IF.
function meansubtracted_IF_Callback(hObject, eventdata, handles)
% hObject    handle to meansubtracted_IF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (isfield(handles,'inst_firing'))
    if (size(handles.cell)>1)
        warndlg('this function not valid for multiple cells.');
        uiwait;
        return;
    end
    window=4;
    motif_length=(handles.motifTemplate(end)-handles.motifTemplate(1)+0.2)*1000; %length in ms
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
    if (isfield(handles,'lineIn'))
        line(handles.xlimit,[handles.lineIn-1 handles.lineIn-1],'color','r');
        line(handles.xlimit,[handles.lineOut handles.lineOut],'color','k');
    end
    if (get(handles.alignment,'Value') == get(handles.alignment,'Max'))
        handles=alignment(handles);
    end
    guidata(hObject, handles);
else
    
    warndlg('First plot the instanteneous firing rate!');
    uiwait;
end

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
    handles.filenums=getAnnotationVector(handles.audioAnnotation,'filenum');
    if isfield(handles,'cell')
        for i=1:length(handles.cell)    
            from(i)=find(handles.filenums>=handles.cell{i}(1,1).filenum & handles.filenums<=handles.cell{i}(1,end).filenum ,1,'first');
        end
        handles.fromIndexAudio=max(from);
        set(handles.from,'String',num2str(handles.filenums(handles.fromIndexAudio)));
        for i=1:length(handles.cell)    
            to(i)=find(handles.filenums>=handles.cell{i}(1,1).filenum & handles.filenums<=handles.cell{i}(1,end).filenum ,1,'last');
        end
        handles.toIndexAudio=min(to);
        set(handles.to,'String',num2str(handles.filenums(handles.toIndexAudio)));
    end
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
subplot(handles.raster);
cla;
subplot(handles.spikeStats)
cla;
if (~isfield(handles,'toIndexAudio') || ~isfield(handles,'fromIndexAudio'))
    warndlg('set the records that you want to display');
    uiwait;
else
    [handles.allSeq handles.startSeq]=sortSequences(handles, handles.numberSyll); %sort the sequences of handles.numSyll in terms of how often they appear

    histx=(1:10^handles.numberSyll); %the histogram of the sequences.
    hist_seq=hist(handles.allSeq,histx);
    %choose the 5 most abundant sequences to display in the sequence pop-up
    %window
    [sortSeq, sortIndex]=sort(hist_seq,'descend');
    numSeq=length(find(sortSeq>0));
    if isempty(numSeq)
        warndlg('no sequences');
        uiwait;
        return;
    elseif (numSeq<5)
        popupNum=numSeq;
    else
        popupNum=7; %display only the four most frequent sequences
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
            popupSeqString(i,:)=[num2str(sortIndex(i)),' (n=', freq(i,:),')'];
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

function [allSeqNumbers startSeq]=sortSequences(handles,numberSyll)
%make a cell array with all the song sequences


[allSeq startSeq counter]=getSequences(handles,numberSyll);
if (allSeq==0)
    warndlg('no sequences in the specified range');
    uiwait;
    return;
end
%go through this cell array and build up an array of sequences
%numberSyll long
for i=1:numberSyll
    multiplier(i)=10^(numberSyll-i);
end
allSeqNumbers=[];
for i=1:counter-1
    allSeqNumbers(i)=dot(allSeq(i,:),multiplier);
end

function [allSeq startSeq counter]=getSequences(handles,numberSyll)
%gets all the sequences containing numSyll
counter=1;

for i=handles.fromIndexAudio:handles.toIndexAudio
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
if ~exist('allSeq')
    allSeq=0;
    startSeq=0;
    counter=0;
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

%fill up the syllable and motif matrices 
subplot(handles.raster);
cla;
subplot(handles.spikeStats)
cla;
for i=1:length(handles.cell)
    handles.spikeIndex{i}=getFieldVector(handles.cell{i},'filenum');%gets the filenumbers that go with the index in the cell files
end
segType=handles.audioAnnotation{1,handles.chosenStartSeq(1,1)}.segType;
handles.sequence=[];
handles.sequence(1:handles.numberSyll)=segType(handles.chosenStartSeq(1,2):handles.chosenStartSeq(1,2)+handles.numberSyll-1);
displayTemplate(handles,handles.sequence);
guidata(hObject, handles);



function displayTemplate(handles,sequence,syllTimes,xdata)

%syllTimes is the average syll duration
if (~isfield(handles,'templates'))
    warndlg('Load a Template');
    uiwait;
else
    noDisp=0;
    fs=handles.audioAnnotation{1,1}.fs;
    templateDisplay(1:floor(fs*0.1))=0; %fill up buffer of 100 ms with zeros
    segTypes=getFieldVector(handles.templates.wavs,'segType');
    for i=1:length(sequence)
            presentSyll=find(segTypes==sequence(i),1,'first');
            if isempty(presentSyll)
               warndlg('no template for the chosen syllable');
               uiwait;
               templateDisplay=0;
               noDisp=1;%don't display anything
            elseif exist('syllTimes')
                noPoints= floor((syllTimes(2*i)-syllTimes(2*i-1))*fs); %number of points in proto syllable
                resampledSyll=resample(handles.templates.wavs(1,presentSyll).wav,noPoints,length(handles.templates.wavs(1,presentSyll).wav));
                templateDisplay(end+1:end+length(resampledSyll))=resampledSyll;
                if (i~=length(sequence))    
                    pause=floor((syllTimes(2*i+1)-syllTimes(2*i))*fs);
                    templateDisplay(end+1:end+pause)=0;
                end

            else    
                templateDisplay(end+1:end+length(handles.templates.wavs(1,presentSyll).wav))=handles.templates.wavs(1,presentSyll).wav;
                templateDisplay(end+1:end+floor(fs*0.015))=0;  
            end
     
    end
    if (noDisp==1)
        axes(handles.specgram);
        cla;
    else
        
        templateDisplay(end+1:end+floor(fs*0.1))=0;
        axes(handles.specgram);
        if (exist('xdata'))    
            displaySpecgramQuick(templateDisplay, fs, [0,10000],[],0,xdata);
        else
            displaySpecgramQuick(templateDisplay, fs, [0,10000],[],0);
        end
        set(handles.specgram,'XTick',[]);
        set(handles.specgram,'YTick',[]);
        xlabel('');
        ylabel('');
    end
end    

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
numberSyllTypes=max(max(allSeq));
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
if isequal(handles.filename,0) || isequal(handles.dirpath,0)
else
    load([handles.dirpath handles.filename]);
    cd (handles.dirpath);
    handles.templates=templates;
    guidata(hObject, handles);
end


% --- Executes on button press in ISI.
function ISI_Callback(hObject, eventdata, handles)
% hObject    handle to ISI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (isfield(handles,'cell') && isfield(handles,'fromIndexAudio'))
    sum_ifr_hist=zeros(201,1);
    if (~isfield(handles,'spikeIndex'))
        handles.spikeIndex=getFieldVector(handles.cell,'filenum');%gets the filenumbers that go with the index in the cell files
    end
    for i=handles.fromIndexAudio:handles.toIndexAudio
         if (length(handles.audioAnnotation{1,i}.segFileStartTimes>2))
             pauses=[];
             for j=1:length(handles.audioAnnotation{1,i}.segFileStartTimes)-1
                 pauses(j)=handles.audioAnnotation{1,i}.segFileStartTimes(j+1)-handles.audioAnnotation{1,i}.segFileEndTimes(j);
                 song_interruptions=find(pauses>0.3); %gives the syllable index after which there is a pause (song interruption, no sound for 300ms)
             end
             presentSpikeIndex=find(handles.spikeIndex{1}==handles.audioAnnotation{1,i}.filenum);
             spikes=handles.cell{1}(1,presentSpikeIndex).spikes;
             if isempty(song_interruptions)
                 spikeIndex = find (spikes>handles.audioAnnotation{1,i}.segFileStartTimes(1)-0.15 & spikes<handles.audioAnnotation{1,i}.segFileEndTimes(end)+0.1);
                 presentSpikes=spikes(spikeIndex);
             else
                 for j=1:length(song_interruptions)
                     if (j==1)
                        spikeIndex = find (spikes>handles.audioAnnotation{1,i}.segFileStartTimes(1)-0.1 & spikes<handles.audioAnnotation{1,i}.segFileEndTimes(song_interruptions(1))+0.1);
                        presentSpikes=spikes(spikeIndex);
                     elseif (j==length(song_interruptions))
                        spikeIndex = find (spikes>handles.audioAnnotation{1,i}.segFileStartTimes(song_interruptions(j)+1)-0.1 & spikes<handles.audioAnnotation{1,i}.segFileEndTimes(end)+0.1);
                        presentSpikes(end+1:end+length(spikeIndex),1)=spikes(spikeIndex);
                     else
                        spikeIndex = find (spikes>handles.audioAnnotation{1,i}.segFileStartTimes(song_interruptions(j)+1)-0.1 & spikes<handles.audioAnnotation{1,i}.segFileEndTimes(song_interruptions(j+1))+0.1); 
                        presentSpikes(end+1:end+length(spikeIndex))=spikes(spikeIndex);
                     end
                 end

             end
             for k=1:length(presentSpikes)
                 presentSpikes(k)=presentSpikes(k)+rand*(1/40000);%rand is there to make the distribution smoother
             end
             if (~isempty(presentSpikes))
                 isi=diff(presentSpikes);                 
                 isiIndex=find(isi<0.2); %get rid of the spikes spanning consecutive motifs
                 isi=isi(isiIndex);
                 num=ones(length(isi),1);
                 if (i==1880)
                     
                     i
                 end
                 ifr=num./isi;
                 ifr_bins=(0:4:800);
                 ifr_hist=zeros(201,1);
                 ifr_hist = histc(ifr,ifr_bins);
                 if (size(ifr_hist)==size(sum_ifr_hist))
                    sum_ifr_hist=sum_ifr_hist+ifr_hist;
                 else
                     warndlg('something is up');
                     uiwait;
                     i
                 end
             end

           
         end         
     end

    sum_ifr_hist = sum_ifr_hist/sum(sum_ifr_hist);% (normalized)
    if (~isfield(handles,'isiFig'))
        handles.isiFig=figure('Name','ISI Distribution','NumberTitle','off');
    else
        figure(handles.isiFig);
    end
    plot(ifr_bins,sum_ifr_hist);
end   
guidata(hObject, handles);




% --- Executes on button press in BeforeAfter.
function BeforeAfter_Callback(hObject, eventdata, handles)
% hObject    handle to BeforeAfter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BeforeAfter
if (length(handles.cell)>1)
    warndlg('this function not valid for multiple cells.');
    uiwait;
    return;
end
if ~isfield(handles,'fromIndexAudio')
    warndlg('Choose records first');
    uiwait;
    return;
end
after=get(hObject,'Value');
chosenSyll = get(handles.syllable, 'Value'); %chosen syllable
[allSeq startSeq]=sortSequences(handles, 2);%gets all the sequences for 2 syllables
spikeIndex=getFieldVector(handles.cell{1},'filenum');
color{1}=[0 1 0];
color{2}=[1 0 0];
color{3}=[0 0 1];
color{4}=[0 0 0];
color{5}=[0 1 1];

if (allSeq==0)
    return;
end
if (after==1)
    chosenSeq=find(allSeq>chosenSyll*10 & allSeq<(chosenSyll+1)*10);%finds the sequences where the chosen sylable if followed by another syllable
else
    allSeq=allSeq-chosenSyll;
    [m,n] = size(allSeq);
    tens(1:m,1:n)=10;
    chosenSeq=find(rem(allSeq,tens)==0);
    allSeq=allSeq+chosenSyll;
end


if isempty(chosenSeq)
    warndlg('no such syllable');
    uiwait;
    return;
else
    chosenSeq=allSeq(chosenSeq);%chosenSeq is an array of all the sequences starting with the chosen syllable
    %choose the tree most prevalent transitions
    histx=(1:100); %the histogram of the sequences.
    hist_seq=hist(chosenSeq,histx);
    [sortSeq, sortIndex]=sort(hist_seq,'descend'); %sortIndex(i) is the ith most numerous sequence sortSeq(i) is the number of instances of the ith sequence
    numSeq=length(find(sortSeq>4));%there has to be at least 2 instances of the syllable transition to be counted here
    if (isempty(numSeq)||numSeq==0)
        warndlg('not enough syllable transitions for the chosen syllable');
        uiwait;
        return;
    elseif (numSeq>5)
        numSeq=5; %display a max of three transitions
    end  
    spikes_missing(1:numSeq)=0;
    counter=0;
    count(1:numSeq)=0;
    for i=1:numSeq
        startIndex=find(allSeq==sortIndex(i)); %finds the index of the chosen sequence as they appear in handles.all_seq
        chosenStartSeq=[];
        for j=1:length(startIndex)
            chosenStartSeq(j,:)=startSeq(startIndex(j),:); 
        end
        for j=1+counter:size(chosenStartSeq,1)+counter

            syll(1:2:3,j)=handles.audioAnnotation{1,chosenStartSeq(j-counter,1)}.segFileStartTimes(chosenStartSeq(j-counter,2):chosenStartSeq(j-counter,2)+1);
            syll(2:2:4,j)=handles.audioAnnotation{1,chosenStartSeq(j-counter,1)}.segFileEndTimes(chosenStartSeq(j-counter,2):chosenStartSeq(j-counter,2)+1);
            presentSpikeIndex=find(spikeIndex==handles.audioAnnotation{1,chosenStartSeq(j-counter,1)}.filenum);
            if isempty(presentSpikeIndex)
                warndlg(['no spikes for file' num2str(handles.audioAnnotation{1,chosenStartSeq(j-counter,1)}.filenum)]);
                spikes{j}=[];
                uiwait;
            else
                spikes{j}=handles.cell{1}(1,presentSpikeIndex).spikes;
                colorIndex(j)=i;
            end
        end
        counter=counter+sortSeq(i);
    end %now we have a 'syll' file with and corresponding 'spikes' sorted by what is followed 

       %display them as spike trains and
        %show their psths.
    noMotifs=size(syll,2);

    index=3; %use this as origin
    syllIndex=1; %display average times
    startBuffer=-0.3;
    finishBuffer=0.1;

    for i=1:noMotifs
        syllTimes(:,i)=syll(:,i)-syll(index,i);
        spikes{i}=spikes{i}-syll(index,i);
    end

    motifTemplate(1)=mean(syllTimes(index,:));
    motifTemplate(2)=mean(syllTimes(index+1,:));
    syll(1,1)=mean(syllTimes(syllIndex,1:sortSeq(1)));
    syll(1,2)=mean(syllTimes(syllIndex+1,1:sortSeq(1)));
    axes(handles.raster);
    cla;
    xlim([startBuffer motifTemplate(2)+finishBuffer]);
    ylim([0,noMotifs]);
    set(gca,'YDir','normal');   
    for j=1:2
        line([syll(1,j) syll(1,j)], [0 sortSeq(1)],'color',[1 0.7 0.7],'LineWidth',0.1);      
    end
    if (numSeq>1)
        for i=2:numSeq
            syll(i,1)=mean(syllTimes(syllIndex,sum(sortSeq(1:i-1))+1:sum(sortSeq(1:i))));
            syll(i,2)=mean(syllTimes(syllIndex+1,sum(sortSeq(1:i-1))+1:sum(sortSeq(1:i))));
            for j=1:2
                line([syll(i,j) syll(i,j)], [sum(sortSeq(1:i-1)) sum(sortSeq(1:i))],'color',[1 0.7 0.7],'LineWidth',0.1);      
            end
        end
    end
 
    for i=1:noMotifs    
        spikeTrain=spikes{i};
%         if (get(handles.timewarping,'Value')==1)
%             markerTimes=syllTimes(1:2:4,i); %warp the syllable
%             spikeTrain = warpSpikeTimes(spikeTrain, markerTimes, [mean(syllTimes(1,:)) mean(syllTimes(3,:)) ]);
%             spikes{i}=spikeTrain; 
%             for j=1:length(spikeTrain)
%                 line([spikeTrain(j) spikeTrain(j)], [i-0.9 i-0.1],'color',color{colorIndex(i)});
%             end
%         else
            for j=1:length(spikeTrain)
            line([spikeTrain(j) spikeTrain(j)], [i-0.9 i-0.1],'color',color{colorIndex(i)});
            end
%         end

    end  
    for i=1:numSeq
        htext(i)=text(motifTemplate(2)+finishBuffer/1.2,sum(sortSeq(1:i))-sortSeq(i)/2,num2str(sortIndex(i)),'HorizontalAlignment','center','Color','k');
    end
    handles.xlimit=get(gca,'xlim');
    handles.ylimit=get(gca,'ylim');
    %make the lines that show where on average the followup/preceeding syllables come
    
    handles.motifTemplate=motifTemplate;
    handles.numberSyll=1;
    line([0 0],handles.ylimit, 'color',[0.2 0 0],'LineWidth',0.1);
     
    for i=1:numSeq   
        line(handles.xlimit,[sum(sortSeq(1:i)) sum(sortSeq(1:i))],'color','r');
    end
     %plot the psth
     subplot(handles.spikeStats);
     cla;
     bin=[];
     bin_size=3; %bin size for psth in ms
     motif_length=(motifTemplate(end)-motifTemplate(1)+0.8)*1000; %length in ms
     for i=1:floor(motif_length/bin_size)
         bin(i)=motifTemplate(1)-0.4+(i*bin_size)/1000;
     end
     bin_size=bin(2)-bin(1);
     psth_seq=zeros(length(bin),numSeq);
     for j=1:numSeq
         if (j==1)
             start=1;
             finish=sortSeq(1);
         else
             start=finish+1;
             finish=start+sortSeq(j)-1;
         end
         for i=start:finish
             spikeTrain=spikes{i};
             if(~isempty(spikeTrain))
                psth = histc(spikeTrain,bin);
                psth_seq(:,j)=psth+psth_seq(:,j);
             else
                 spikes_missing(j)=spikes_missing(j)+1;
             end
         end
         psth_seq(:,j)=psth_seq(:,j)/((sortSeq(j)-spikes_missing(j))*bin_size);
         plot(bin,psth_seq(:,j),'color',color{j});
         hold on;
         xlim([(motifTemplate(1)+startBuffer) (motifTemplate(end)+finishBuffer)]);

     end

end
if (after==1)
    xdata=[mean(syll(:,1))-0.1 mean(syll(:,2))+0.1];
else
    xdata=[motifTemplate(1)-0.1 motifTemplate(end)+0.1];
end
displayTemplate(handles,chosenSyll,motifTemplate,xdata);
xlim([(motifTemplate(1)+startBuffer) (motifTemplate(end)+finishBuffer)]);



    
% --- Executes on selection change in syllable.
function syllable_Callback(hObject, eventdata, handles)
% hObject    handle to syllable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns syllable contents as cell array
%        contents{get(hObject,'Value')} returns selected item from syllable


% --- Executes during object creation, after setting all properties.
function syllable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to syllable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


