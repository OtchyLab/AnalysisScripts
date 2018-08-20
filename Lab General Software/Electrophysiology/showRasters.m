 
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

% Last Modified by GUIDE v2.5 22-Sep-2010 15:03:40

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
set(handles.raster, 'ButtonDownFcn', @cb_raster_click);
set(handles.specgram,'XTick',[]);
set(handles.specgram,'YTick',[]);
set(handles.raster,'XTick',[]);
%set(handles.raster,'YTick',[]);
set(handles.popupSequence,'Value',1);
handles.fs=44100; %sampling rate
handles.eventLine1=[];
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
handles.treshold=80;%default value for the treshold I don't really know where to put it so it is here
[handles.filename,handles.dirpath]=uigetfile('*.mat');
handles.name_cell=handles.filename(1:end-4);
if isequal(handles.filename,0) || isequal(handles.dirpath,0)
else
    if isfield(handles,'cell') %if there already exists a cell, add to it or replace it?
        button = questdlg('Do you want to add this cells to previously loaded cells (simultaneously recorded)?','Add a cell','Add','Replace','Add');
        if (strcmp(button,'Replace'))
            handles.cell={};
            handles.psth_drug=[];
            handles.psth_nodrug=[];
            handles.drugFiles=[];
            handles.isiDrugCells=[];
            handles.isiNoDrugCells=[];
            handles.isiCellsDrugTimeWeighted=[];
            handles.isiCellsNoDrugTimeWeighted=[];
            handles.activeDrug=[];
            handles.activeNoDrug=[];
            handles.frDrug=[];
            handles.frNoDrug=[];

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
if isfield(handles,'audioAnnot')
        
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
        handles.audioAnnot = audioAnnotation.elements;
        handles.audioKeys = audioAnnotation.keys;
        handles.filenums=getAnnotationVector(handles.audioAnnot,'filenum');
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
    color{1}=[0 0 1];
    color{2}=[0 0 0];
    color{3}=[1 0 0];
    color{4}=[0 1 0];
    color{5}=[0 1 1];
    [handles.spikes handles.motif handles.syll]=fillUpRaster(handles);%gives back a spike array, and syllable times
    
    %For data extraction:
    keys = {};
    times = handles.motif;
    range = [times(1,:)-.25;times(2,:)+.25];
    
    %This section is to excerpt the data for Haim's analysis
    for i=1:length(handles.chosenStartSeq)
        keys{i} = handles.keys{handles.chosenStartSeq(i,1)};
        [HWChannels, tempaudio, time, startSamp] = daq_readDatafileBence(keys{i});

        audio{i} = tempaudio(round(range(1,i)*44150):round(range(2,i)*44150));

        tempspikes = handles.spikes{i};
        spikes{i} = tempspikes(tempspikes>=range(1,i) & tempspikes<=range(2,i))-range(1,i);
    end
    
    
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
            from=(i-1)*length(handles.cell)+(j)-0.9;
            to=(i-1)*length(handles.cell)+j-0.1;
            spikeIndex=find(handles.spikes{j,i}>=(handles.motifTemplate(1)-0.05) & handles.spikes{j,i}<=(handles.motifTemplate(end)));
            spikes=handles.spikes{j,i}(spikeIndex);
              
            temp=find(spikes<0.4);
            handles.spikeCount{j,i}=length(temp);
            
            if (get(handles.spikewarping,'Value')==0)
                for k=1:length(spikes)
                    line([spikes(k) spikes(k)], [from to], 'color', color{j});
                end
            end
        end
    end    
    %calculate the FR in each rendition
%      [drugOut drugIn drug nodrug]=getDrugInfo(handles);
%     for i=1:length(nodrug)
%         noDrugFR(i)=handles.spikeCount{1,nodrug(i)}/handles.motifTemplate(end)-(handles.motifTemplate(1)-0.05);
%     end
%     for i=1:length(drug)
%         DrugFR(i)=handles.spikeCount{1,drug(i)}/handles.motifTemplate(end)-(handles.motifTemplate(1)-0.05);
%     end
%     if length(drug)>0 & length(nodrug)>0
%         lengthDrug=length(drug)
%         lengthNoDrug=length(nodrug)
%         mND=mean(noDrugFR)
%         mD=mean(DrugFR)
%         sND=std(noDrugFR)
%         sD=std(DrugFR)
%         [h,p] = ttest2(noDrugFR,DrugFR)
%     end
    handles.xlimit=get(gca,'xlim');
    handles.ylimit=get(gca,'ylim');
    if (get(handles.alignment,'Value') == get(handles.alignment,'Max'))
        handles=alignment(handles);
    end
    
    %plot the psth
     handles.bin=[];
     bin_size=3; %bin size for psth in ms
     motif_length=(handles.motifTemplate(end)-handles.motifTemplate(1)+0.04)*1000; %length in ms
     for i=1:floor(motif_length/bin_size)
         handles.bin(i)=handles.motifTemplate(1)-0.04+(i*bin_size)/1000;
     end
     bin_size=handles.bin(2)-handles.bin(1);
     if (get(handles.spikewarping,'Value')==1)
        cellRef=1;%uses cell 1 as reference
        spikes=cat(1,handles.spikes{cellRef,1:handles.noMotifs}); 
        psth = histc(spikes,handles.bin)/(handles.noMotifs*bin_size);
        startTime=handles.motifTemplate(1)-0.04;
        endTime=handles.motifTemplate(end);
        sampleRate=handles.fs;
        GaussWidth=0.008; %in seconds
        motifShifts = crossCorrelate(1,handles.noMotifs, handles.spikes,sampleRate, 0, 0, GaussWidth, startTime,endTime,floor(1/bin_size),psth, cellRef);

        for j=1:length(handles.cell)
            for i=1:handles.noMotifs 
                handles.spikes{j,i}=handles.spikes{j,i}+ motifShifts(i)/1000;
            end
        end
        % shift the spike train according to the lag that gives the
        % highest corr.
        subplot(handles.raster);
        cla;
        xlim([(handles.motifTemplate(1)-0.1) (handles.motifTemplate(end)+0.1)]);
        ylim([0,handles.noMotifs]);
        set(gca,'YDir','normal');
        for j=1:length(handles.cell)
            for i=1:handles.noMotifs
                from=(i-1)*length(handles.cell)+(j)-0.9;
                to=(i-1)*length(handles.cell)+j-0.1;
                spikeIndex=find(handles.spikes{j,i}>=(handles.motifTemplate(1)-0.1) & handles.spikes{j,i}<=(handles.motifTemplate(end)+0.1));
                spikes=handles.spikes{j,i}(spikeIndex);
              
                    
          
                for k=1:length(spikes)
                    line([spikes(k) spikes(k)], [from to], 'color', color{j});
                end
            end
        end        
     end
     % check wether there has been a drug infusion
     if (length(handles.cell)==1)
         [drugOut drugIn drug nodrug]=getDrugInfo(handles);
         if (~isempty(drugIn))
             if (length(drugIn)>length(drugOut))
                 drugOut(end+1)=max(handles.chosenStartSeq(:,1));
             end        
             for i=1:length(drugIn)
                handles.lineIn(i)=find(handles.chosenStartSeq(:,1)>=drugIn(i),1,'first');
                if isempty(handles.lineIn)
                    return;
                end
                subplot(handles.raster);
                line(handles.xlimit,[handles.lineIn(i)-1 handles.lineIn(i)-1],'color','r');
             end
         end
         if(~isempty(drugOut))
             for i=1:length(drugOut)
                    if (~isempty(find(handles.chosenStartSeq(:,1)<=drugOut(i),1,'last')));
                        handles.lineOut(i)=find(handles.chosenStartSeq(:,1)<=drugOut(i),1,'last');
                        line(handles.xlimit,[handles.lineOut(i) handles.lineOut(i)],'color','k');   
                    end
             end
         end
         handles.bin=[];
         bin_size=3; %bin size for psth in ms
         motif_length=(handles.motifTemplate(end)-handles.motifTemplate(1)+0.8)*1000; %length in ms
         for i=1:floor(motif_length/bin_size)
             handles.bin(i)=handles.motifTemplate(1)-0.4+(i*bin_size)/1000;
         end
         bin_size=handles.bin(2)-handles.bin(1);
         handles.psth_drug=zeros(1,length(handles.bin));
         handles.psth_nodrug=zeros(1,length(handles.bin));
         [cND, handles.nodrugIndex] = ismember((handles.chosenStartSeq(:,1)), nodrug','rows');
         ND=find(cND==1);
         handles.nodrugIndex=ND;
         nodrugSpikes=cat(1,handles.spikes{1,handles.nodrugIndex});
         [cD, handles.drugIndex] = ismember((handles.chosenStartSeq(:,1)), drug','rows');
         D=find(cD==1);
         handles.drugIndex=D;
         drugSpikes=cat(1,handles.spikes{1,handles.drugIndex});
         handles.psth_drug = hist(drugSpikes,handles.bin)/(length(D)*bin_size);
         handles.psth_nodrug = hist(nodrugSpikes,handles.bin)/(length(ND)*bin_size);
         %calculate the number of spikes for the 
         
         first=find(handles.bin>(handles.motifTemplate(1)-0.05),1,'first');
         last=find(handles.bin<(handles.motifTemplate(end)),1,'last');  
         if ~isempty(handles.psth_drug)
            handles.frDrug=sum(handles.psth_drug(first:last))/(last-first);
         end
         if ~isempty(handles.psth_nodrug)
            handles.frNoDrug=sum(handles.psth_nodrug(first:last))/(last-first);
         end
         
         
         if length(ND)>8
            for i=1:length(ND)
                spikesND{i}=handles.spikes{1,handles.nodrugIndex(i)}(find(handles.spikes{1,handles.nodrugIndex(i)}>(handles.motifTemplate(1)-0.05) & (handles.spikes{1,handles.nodrugIndex(i)}<(handles.motifTemplate(end)))));
                handles.numSpikesND(i)=length(spikesND{i});
                ISI=diff(spikesND{i});
                ISI_HF=find(ISI<0.0025001);%index of spikes that are in bursts
                spikes_adj=length(find(diff(ISI_HF)==1)); 
                HF_spikes=2*length(ISI_HF)-spikes_adj;
                fractionND(i)=HF_spikes/length(spikesND{i});
                FRnodrug(i)=handles.numSpikesND(i)/(handles.motifTemplate(end)-(handles.motifTemplate(1)-0.05));
            end    
         end
         FRND=mean(FRnodrug);
         if length(D)>8
            for i=1:length(D)
                spikesD{i}=handles.spikes{1,handles.drugIndex(i)}(find(handles.spikes{1,handles.drugIndex(i)}>(handles.motifTemplate(1)-0.05) & (handles.spikes{1,handles.drugIndex(i)}<(handles.motifTemplate(end)))));
                handles.numSpikesD(i)=length(spikesD{i});
                ISI=diff(spikesD{i});
                ISI_HF=find(ISI<0.0025001);%index of spikes that are in bursts
                spikes_adj=length(find(diff(ISI_HF)==1)); 
                HF_spikes=2*length(ISI_HF)-spikes_adj;
                fractionD(i)=HF_spikes/length(spikesD{i});
                FRdrug(i)=handles.numSpikesD(i)/(handles.motifTemplate(end)-(handles.motifTemplate(1)-0.05));
            end 
         end
         %is the firing rate significantly different for the chosen song motif?
         if length(D)>8 && length(ND)>8   
            [handles.sign,handles.pFR] = ttest2(FRnodrug,FRdrug);
            [sign,p]=ttest2(FRnodrug,FRdrug);
            [signFraction, pfraction]=ttest2(fractionND,fractionD)
            FRD=mean(FRdrug);
            FRND=mean(FRnodrug);
            FD=mean(fractionD)
            FND=mean(fractionND)
            FDSD=std(fractionD)
            FNDSD=std(fractionND)
            SDD=std(FRdrug);
            SDND=std(FRnodrug);
            
            NumberDrug=length(D);
            NumberNodrug=length(ND);
            c=corrcoef(handles.psth_nodrug(first:last),handles.psth_drug(first:last));
            handles.PSTHcorr=c(1,2);
         else
            handles.PSTHcorr=[];
            handles.numSpikesD=[];
            handles.numSpikesND=[];
         end
         %add a function that calculates sparseness during no drug
         %condition
         sumFR=sum(handles.psth_nodrug(first:last));
         BigN=length(handles.psth_nodrug(first:last));
        
         handles.motifFRND=FRND;
         for i=1:BigN
            p(i)=handles.psth_nodrug(first+i-1)/sumFR;
            if p(i)>0
                sparseness(i)=(p(i)*log(p(i)));
            else
                sparseness(i)=0;
            end
         end
         handles.sparseness=1+(sum(sparseness)/log(BigN));
     
         
         subplot(handles.spikeStats);
         cla;
         plot(handles.bin,handles.psth_nodrug,handles.bin,handles.psth_drug,'r');
         xlim([(handles.motifTemplate(1)-0.1) (handles.motifTemplate(end)+0.1)]);
         
         %do it also for bin_size 10 ms
         bin_size=10;
         handles.bin=[];
         handles.psth_nodrug=[];
         for i=1:floor(motif_length/bin_size)
             handles.bin(i)=handles.motifTemplate(1)-0.4+(i*bin_size)/1000;
         end
         bin_size=handles.bin(2)-handles.bin(1);
         handles.psth_nodrug=zeros(1,length(handles.bin));
         [cND, handles.nodrugIndex] = ismember((handles.chosenStartSeq(:,1)), nodrug','rows');
         ND=find(cND==1);
         handles.nodrugIndex=ND;
         nodrugSpikes=cat(1,handles.spikes{1,handles.nodrugIndex});
         handles.psth_nodrug = hist(nodrugSpikes,handles.bin)/(length(ND)*bin_size);
         first=find(handles.bin>(handles.motifTemplate(1)-0.05),1,'first');
         last=find(handles.bin<(handles.motifTemplate(end)),1,'last');  
          sumFR=sum(handles.psth_nodrug(first:last));
         BigN=length(handles.psth_nodrug(first:last));
         p=[];
         
         for i=1:BigN
            p(i)=handles.psth_nodrug(first+i-1)/sumFR;
            if p(i)>0
                sparseness10(i)=(p(i)*log(p(i)));
            else
                sparseness10(i)=0;
            end
         end
         handles.sparseness10=1+(sum(sparseness10)/log(BigN));
       
     else
         %if more than 1 cell
         motifs=1:handles.noMotifs;
         handles.bin=[];
         bin_size=3; %bin size for psth in ms
         motif_length=(handles.motifTemplate(end)-handles.motifTemplate(1)+0.8)*1000; %length in ms
         for i=1:floor(motif_length/bin_size)
             handles.bin(i)=handles.motifTemplate(1)-0.4+(i*bin_size)/1000;
         end
         bin_size=handles.bin(2)-handles.bin(1);
         subplot(handles.spikeStats);
         cla;
         for j=1:length(handles.cell)
            spikes=cat(1,handles.spikes{j,motifs});
            psth=histc(spikes,handles.bin)/(handles.noMotifs*bin_size);
            plot(handles.bin,psth,'color',color{j});
            xlim([(handles.motifTemplate(1)-0.1) (handles.motifTemplate(end)+0.1)]);
            hold on;
         end
         
     end
end

guidata(hObject, handles);

function [drugOut drugIn drug nodrug]=getDrugInfo(handles)
 drugStatus=getAnnotationVector(handles.audioAnnot,'drugindex');
 diffDrugStatus=diff(drugStatus);
 drugIn=find(diffDrugStatus>0)+1; %gives the filenumber at which the switch occurs 
 drug=find(drugStatus~=1);
 drug=drug(find(drug>=handles.fromIndexAudio & drug<=handles.toIndexAudio));
 nodrug=find(drugStatus==1);
 nodrug=nodrug(find(nodrug>=handles.fromIndexAudio & nodrug<=handles.toIndexAudio));
 drugInIndex=find(drugIn>=handles.fromIndexAudio & drugIn<=handles.toIndexAudio);
 drugIn=drugIn(drugInIndex);
 drugOut=find(diffDrugStatus<0);
 drugOutIndex=find(drugOut>=handles.fromIndexAudio & drugOut<=handles.toIndexAudio);
 drugOut=drugOut(drugOutIndex);
         

function [spikes motif syll]=fillUpRaster(handles)

for i=1:size(handles.chosenStartSeq,1)
    syll(1:2:handles.numberSyll*2-1,i)=handles.audioAnnot{1,handles.chosenStartSeq(i,1)}.segFileStartTimes(handles.chosenStartSeq(i,2):handles.chosenStartSeq(i,2)+handles.numberSyll-1);
    syll(2:2:handles.numberSyll*2,i)=handles.audioAnnot{1,handles.chosenStartSeq(i,1)}.segFileEndTimes(handles.chosenStartSeq(i,2):handles.chosenStartSeq(i,2)+handles.numberSyll-1);
    motif(1,i)= handles.audioAnnot{1,handles.chosenStartSeq(i,1)}.segFileStartTimes(handles.chosenStartSeq(i,2));
    motif(2,i)= handles.audioAnnot{1,handles.chosenStartSeq(i,1)}.segFileEndTimes(handles.chosenStartSeq(i,2)+handles.numberSyll-1);
    for j=1:length(handles.cell)
        presentSpikeIndex=find(handles.spikeIndex{j}==handles.audioAnnot{1,handles.chosenStartSeq(i,1)}.filenum);  
        if isempty(presentSpikeIndex)
            warndlg(['no spikes for file ' num2str(handles.audioAnnot{1,handles.chosenStartSeq(i,1)}.filenum) ' cell number' num2str(j)]);
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
    for cellNo=1:length(handles.cell)
%         if (length(handles.cell)>1)
%             warndlg('this function not valid for multiple cells.');
%             uiwait;
%             return;
%         end
        startTime=handles.motifTemplate(1)-0.03;
        endTime=handles.motifTemplate(end)+0.02;
        sampleRate=handles.fs;
        GaussWidth=0.008; %in seconds
        for i=1:size(handles.spikes,2)
            spikes{i}=handles.spikes{cellNo,i};
        end
        handles.corr = pairWiseCorrelations(1,handles.noMotifs, spikes,sampleRate, 0, 0, GaussWidth, startTime,endTime);
        %handles.corr=handles.corr(handles.fromIndexAudio:handles.to_index,handles.fromIndexAudio:handles.to_index);
        if (length(handles.cell)==1)
            subplot(handles.raster);
            cla;
            hold off;
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
            handles.correlation=sum_corr/((m*n)-m)
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
            sum_corr=sum(sum(handles.corr))-sum(diag(handles.corr));
            [m,n] = size(handles.corr);
            handles.correlation=sum_corr/((m*n)-m)
        end
    end
else
    warndlg('plot the spikes first')
    uiwait;
end

% tags=input('unique date and cell id','s');
% savedir = 'C:\Users\Tim\Desktop\Data';
% filename=[savedir '\Pur534_data_' tags '.mat'];
% save(filename,'handles');




% --- Executes on button press in inst_firing.
function inst_firing_Callback(hObject, eventdata, handles)
% hObject    handle to inst_firing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (length(handles.cell)>1)
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

if (length(handles.cell)>1)
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
sum_victor=sum(sum(victor))-sum(diag(victor));
victor_sum=sum_victor/((m*n)-m)
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
    if (length(handles.cell)>1)
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
    %The follow ing line was put in 11/5/11 to deal with lotsa of
    %unlabelled syllables; turn this toggle on/off as needed.
    handles.audioAnnot = audioAnnotation.elements;
    %[handles.audioAnnot] = unlabelsylremap(audioAnnotation.elements);
    
    handles.keys = audioAnnotation.keys;
    handles.filenums=getAnnotationVector(handles.audioAnnot,'filenum');
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
    elseif (numSeq<10)
        popupNum=numSeq;
    else
        popupNum=9; %display only the four most frequent sequences
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
song_interruptions=[];

for i=handles.fromIndexAudio:handles.toIndexAudio
    for j=1:length(handles.audioAnnot{1,i}.segFileStartTimes)-1
        pauses(j)=handles.audioAnnot{1,i}.segFileStartTimes(j+1)-handles.audioAnnot{1,i}.segFileEndTimes(j);
        song_interruptions=find(pauses>0.1); %gives the syllable index after which there is a pause (song interruption, no sound for 100ms)

    end

    
    if (length(handles.audioAnnot{1,i}.segFileStartTimes)>=numberSyll)
        for j=1:length(handles.audioAnnot{1,i}.segFileStartTimes)-numberSyll+1
            current_seq(1:numberSyll)=handles.audioAnnot{1,i}.segType(j:j+numberSyll-1);
            if (isempty(find(current_seq<1 | current_seq>9)) && isempty(find(song_interruptions>j-1 & song_interruptions<j+numberSyll-1)))
            %if (isempty(find(current_seq<-1 | current_seq>9)) && isempty(find(song_interruptions>j-1 & song_interruptions<j+numberSyll-1)))
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
handles.Seq=chosenSeq;
startIndex=find(handles.allSeq==chosenSeq); %finds the index of the chosen sequence as they appear in handles.all_seq
handles.chosenStartSeq=[];%gives the annotation number and the number f syllable within that annotation
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
segType=handles.audioAnnot{1,handles.chosenStartSeq(1,1)}.segType;
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
    fs=handles.audioAnnot{1,1}.fs;
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
    
    [drugOut drugIn drug nodrug]=getDrugInfo(handles);
    firingRateNoDrug=[];
    firingRateND=[];
    firingRateDrug=[];
    firingRate=[];
    frTime=[];
    time=0;
    handles.spiketrainsNoDrugs=[];
    handles.spiketrainsDrugs=[];
    for cellNo = 1:length(handles.cell)
        sum_ifr_nodrugHist=zeros(200,1);
        sum_ifr_drugHist=zeros(200,1);
        if (~isfield(handles,'spikeIndex'))
            handles.spikeIndex{cellNo}=getFieldVector(handles.cell{cellNo},'filenum');%gets the filenumbers that go with the index in the cell files
        end
        cellFilenum=getAnnotationVector(handles.audioAnnot,'filenum');
        handles.drugFiles=cellFilenum(drug);
        songbout_noDrugs=1;
        songbout_Drugs=1;
        counter=0;
        
        for i=handles.fromIndexAudio:handles.toIndexAudio %loop over audio records
             spikes=[];
             presentSpikeIndex=[];
             presentSpikes=[];
            if (length(handles.audioAnnot{1,i}.segFileStartTimes>2))
                 pauses=[];
                 for j=1:length(handles.audioAnnot{1,i}.segFileStartTimes)-1
                     pauses(j)=handles.audioAnnot{1,i}.segFileStartTimes(j+1)-handles.audioAnnot{1,i}.segFileEndTimes(j);
                     song_interruptions=find(pauses>0.2); %gives the syllable index after which there is a pause (song interruption, no sound for 200ms)
                 end
                 
                 if(length(handles.audioAnnot{1,i}.segFileStartTimes)==1)
                     song_interruptions=[];
                 end
                 presentSpikeIndex=find(handles.spikeIndex{cellNo}==handles.audioAnnot{1,i}.filenum);
                 if (~isempty(presentSpikeIndex))
                    spikes(1:length(handles.cell{cellNo}(1,presentSpikeIndex).spikes),1)=handles.cell{cellNo}(1,presentSpikeIndex).spikes;
                 end
                 
                 if isempty(song_interruptions)
                     spikeIndex = find (spikes>handles.audioAnnot{1,i}.segFileStartTimes(1)-0.05 & spikes<handles.audioAnnot{1,i}.segFileEndTimes(end));
                     presentSpikes=spikes(spikeIndex);
                     if (isempty(find(handles.drugFiles==handles.audioAnnot{1,i}.filenum)) && (length(handles.cell)==1)) %if no drugs
                        handles.spiketrainsNoDrugs{songbout_noDrugs}=presentSpikes;
                        songbout_noDrugs=songbout_noDrugs+1;
                        firingRateND (end+1)=length(presentSpikes)/(handles.audioAnnot{1,i}.segFileEndTimes(end)-handles.audioAnnot{1,i}.segFileStartTimes(1)+0.05);
                     else
                        handles.spiketrainsDrugs{songbout_Drugs}=presentSpikes;
                        songbout_Drugs=songbout_Drugs+1;
                     end
                     firingRate(end+1)=length(presentSpikes)/(handles.audioAnnot{1,i}.segFileEndTimes(end)-handles.audioAnnot{1,i}.segFileStartTimes(1)+0.05);
                     frTime(end+1)=time+handles.audioAnnot{1,i}.segFileEndTimes(end)-handles.audioAnnot{1,i}.segFileStartTimes(1)+0.05;
                 else
                     for j=1:length(song_interruptions)
                         if (j==1)
                            spikeIndex = find (spikes>handles.audioAnnot{1,i}.segFileStartTimes(1)-0.05 & spikes<handles.audioAnnot{1,i}.segFileEndTimes(song_interruptions(1)));
                            presentSpikes=spikes(spikeIndex);
                             if (isempty(find(handles.drugFiles==handles.audioAnnot{1,i}.filenum)) && (length(handles.cell)==1)) %if no drugs
                                handles.spiketrainsNoDrugs{songbout_noDrugs}=presentSpikes;
                                songbout_noDrugs=songbout_noDrugs+1;
                             else
                                handles.spiketrainsDrugs{songbout_Drugs}=presentSpikes;
                                songbout_Drugs=songbout_Drugs+1;
                             end
                            time=handles.audioAnnot{1,i}.segFileEndTimes(song_interruptions(1))-handles.audioAnnot{1,i}.segFileStartTimes(1)+0.05;
                         elseif (j==length(song_interruptions))
                            spikeIndex = find (spikes>handles.audioAnnot{1,i}.segFileStartTimes(song_interruptions(j)+1)-0.05 & spikes<handles.audioAnnot{1,i}.segFileEndTimes(end));
                            presentSpikes(end+1:end+length(spikeIndex),1)=spikes(spikeIndex);
                            if (isempty(find(handles.drugFiles==handles.audioAnnot{1,i}.filenum)) && (length(handles.cell)==1)) %if no drugs
                                handles.spiketrainsNoDrugs{songbout_noDrugs}=presentSpikes;
                                songbout_noDrugs=songbout_noDrugs+1;
                            else
                                handles.spiketrainsDrugs{songbout_Drugs}=presentSpikes;
                                songbout_Drugs=songbout_Drugs+1;
                            end
                            time=time+handles.audioAnnot{1,i}.segFileEndTimes(end)-handles.audioAnnot{1,i}.segFileStartTimes(song_interruptions(j)+1)+0.05;
                         else
                            spikeIndex = find (spikes>handles.audioAnnot{1,i}.segFileStartTimes(song_interruptions(j)+1)-0.05 & spikes<handles.audioAnnot{1,i}.segFileEndTimes(song_interruptions(j+1))); 
                            presentSpikes(end+1:end+length(spikeIndex),1)=spikes(spikeIndex);
                            if (isempty(find(handles.drugFiles==handles.audioAnnot{1,i}.filenum)) && (length(handles.cell)==1)) %if no drugs
                                handles.spiketrainsNoDrugs{songbout_noDrugs}=presentSpikes;
                                songbout_noDrugs=songbout_noDrugs+1;
                            else
                                handles.spiketrainsDrugs{songbout_Drugs}=presentSpikes;
                                songbout_Drugs=songbout_Drugs+1;
                            end
                            time=time+handles.audioAnnot{1,i}.segFileEndTimes(song_interruptions(j+1))-handles.audioAnnot{1,i}.segFileStartTimes(song_interruptions(j)+1)+0.05;
                         end
                     end
                    firingRate(end+1)=(length(presentSpikes)/time);
                    frTime(end+1)=time;
                 end
                 
                 for k=1:length(presentSpikes)
                     presentSpikes(k,1)=presentSpikes(k,1)+rand*(1/handles.fs);%rand is there to make the distribution smoother
                     
                 end
                 if (~isempty(presentSpikes))
                     isi(1:length(presentSpikes)-1,1)=diff(presentSpikes);                 
                     isiIndex= isi<0.1; %get rid of the spikes spanning consecutive motifs
                     isi=isi(isiIndex);
                     num=ones(length(isi),1);
                     ifr=num./isi;
                     handles.ifrBins=(2:4:798);
                     ifr_hist=zeros(200,1);
                     ifr_hist(:,1) = histc(ifr,handles.ifrBins);
                     spiketimes=round(((presentSpikes-min(presentSpikes))*1000));
                     timeseries=zeros(1,max(max(spiketimes),1001));
                     timeseries(1,spiketimes(2:end))=1;
                     [ACF,lags]=xcorr(timeseries,1000);
                     ACF(1001)=0;
                     if ~exist('sumSize')
                        sumACF=size(presentSpikes,1)*ACF;
                        sumSize=size(presentSpikes,1);
                     else
                         sumACF=sumACF+size(presentSpikes,1)*ACF;
                         sumSize=sumSize+size(presentSpikes,1);
                     end
                     
                     if (i==handles.toIndexAudio)
                         handles.ACF=sumACF/sumSize;
                     end
         
                     if (isempty(find(handles.drugFiles==handles.audioAnnot{1,i}.filenum)) && (length(handles.cell)==1))
                        sum_ifr_nodrugHist=sum_ifr_nodrugHist+ifr_hist;
                        firingRateNoDrug(end+1:end+1+floor(frTime(end)))=firingRate(end); %this makes no sense
                        
                        firingRateNoDrug(end+1)=firingRate(end);
                     else
                        sum_ifr_drugHist=sum_ifr_drugHist+ifr_hist; %if multiple cells do not distinguish between drug and no drug
                        firingRateDrug(end+1)= firingRate(end);
                        firingRateDrug(end+1:end+1+floor(frTime(end)))=firingRate(end);
                        floor(frTime(end))
                        firingRate(end)
                     end

                 end


             end         
        end
        firingRateNDmean=mean(firingRateND)
       
%         term1=(handles.ifrBins * sum_ifr_drugHist)/sum(sum_ifr_drugHist);
%         term1=term1^2;
%          for i=1:length(handles.ifrBins)
%             temp(i)=(handles.ifrBins(i)^2*sum_ifr_drugHist(i));
%          end
%         term2=sum(temp)/sum(sum_ifr_drugHist);
%           
%         activityIndexDrug=term1/term2;
        sum_ifr_drugHist = sum_ifr_drugHist/sum(sum_ifr_drugHist);% (normalized)
        sum_ifr_nodrugHist = sum_ifr_nodrugHist/sum(sum_ifr_nodrugHist);
        handles.isiDrugCells(cellNo,1:200)=sum_ifr_drugHist;
        handles.isiNoDrugCells(cellNo,1:200)=sum_ifr_nodrugHist;
        time(1:200,1)=1./handles.ifrBins';
        handles.isiCellsDrugTimeWeighted(cellNo,1:200)=sum_ifr_drugHist.*time;
        handles.isiCellsNoDrugTimeWeighted(cellNo,1:200)=sum_ifr_nodrugHist.*time;
        for k=1:200
            handles.activeDrug(cellNo,k)=sum(handles.isiCellsDrugTimeWeighted(cellNo,k:200))/sum(handles.isiCellsDrugTimeWeighted(cellNo,1:200));
            handles.activeNoDrug(cellNo,k)=sum(handles.isiCellsNoDrugTimeWeighted(cellNo,k:200))/sum(handles.isiCellsNoDrugTimeWeighted(cellNo,1:200));
        end
    end
    handles.frTimeBins=(1:2:99);
%     if (~isempty(firingRateDrug))
% %         handles.frDrug(cellNo,1:50)=histc(firingRateDrug,handles.frTimeBins);
%     else
% %         handles.frDrug(cellNo,1:50)=0;
%     end
%     if (~isempty(firingRateNoDrug))
% %         handles.frNoDrug(cellNo,1:50)=histc(firingRateNoDrug,handles.frTimeBins);
%     else
% %         handles.frNoDrug(cellNo,1:50)=0;
%     end
%     if (isempty(eventdata))
%         color{1}=[0 0 1];
%         color{2}=[0 0 0];
%         color{3}=[1 0 0];
%         color{4}=[0 1 0];
%         color{5}=[0 1 1];
%         if (~isfield(handles,'isiFig') || ~isfield(handles,'activeFig') || ~isfield(handles,'FRFig'))
%            handles.isiFig=figure('Name','ISI Distribution','NumberTitle','off');  
%            handles.activeFig=figure('Name','%active','NumberTitle','off');
%            handles.FRFig=figure('Name','Firing Rate','NumberTitle','off'); 
%         end

%         if (length(handles.cell)>1)  
%             for cellNo=1:length(handles.cell)
%                 figure(handles.isiFig);
%                 plot(handles.ifrBins,handles.isiDrugCells(cellNo,1:200),'color',color{cellNo});
%                 hold on;
%                 %semilogy(handles.ifrBins,handles.isiCellsTimeWeighted(cellNo,1:200),'color',color{cellNo});
%                 figure(handles.activeFig);
%                 semilogy(handles.ifrBins,handles.activeDrug(cellNo,1:200),'color',color{cellNo});
%                 hold on;
%                 figure(handles.FRFig);
%                 plot(handles.frTimeBins,handles.frDrug(cellNo,1:50),'color',color{cellNo});
%                 hold on;
% 
%             end
% 
%         else
%                 figure(handles.isiFig);
%                 plot(handles.ifrBins,handles.isiDrugCells(1,1:200),'color','r');
%                 hold on;
%                 plot(handles.ifrBins,handles.isiNoDrugCells(1,1:200));
%                 figure(handles.activeFig);
%                 semilogy(handles.ifrBins,handles.activeDrug(1,1:200),'color','r');
%                 hold on;
%                 semilogy(handles.ifrBins,handles.activeNoDrug(1,1:200));
%                 figure(handles.FRFig);
%                 plot(handles.frTimeBins,handles.frDrug(cellNo,1:50),'color','r');
%                 hold on;
%                 plot(handles.frTimeBins,handles.frNoDrug(cellNo,1:50));
%         end
% 
%         hold off;
%         figure(handles.isiFig);
%         hold off;
%         figure(handles.FRFig);
%         hold off;
%      end
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

            syll(1:2:3,j)=handles.audioAnnot{1,chosenStartSeq(j-counter,1)}.segFileStartTimes(chosenStartSeq(j-counter,2):chosenStartSeq(j-counter,2)+1);
            syll(2:2:4,j)=handles.audioAnnot{1,chosenStartSeq(j-counter,1)}.segFileEndTimes(chosenStartSeq(j-counter,2):chosenStartSeq(j-counter,2)+1);
            presentSpikeIndex=find(spikeIndex==handles.audioAnnot{1,chosenStartSeq(j-counter,1)}.filenum);
            if isempty(presentSpikeIndex)
                warndlg(['no spikes for file' num2str(handles.audioAnnot{1,chosenStartSeq(j-counter,1)}.filenum)]);
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
     bin_size=3;%bin size for psth in ms
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




% --- Executes on button press in spikewarping.
function spikewarping_Callback(hObject, eventdata, handles)
% hObject    handle to spikewarping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in cellFile.
function cellFile_Callback(hObject, eventdata, handles)
% hObject    handle to cellFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (length(handles.cell)>1)
    warndlg('this function not valid for multiple cells.');
    uiwait;
    return;
end
% spontSpikesNoDrug=[];
% spontSpikesDrug=[];

% cellFile.spontSpikesNoDrug=[];
% cellFile.spontSpikesDrug=[];
% cellFile.Seq=handles.Seq;
% if ~(isfield(handles,'lineIn'))
%     cellFile.DrugIn=[];
% else
%     cellFile.DrugIn=handles.lineIn;
% end
% if ~(isfield(handles,'lineOut'))
%     cellFile.DrugOut=[];
% else
%     cellFile.DrugOut=handles.lineOut;
% end

%prompt = {'Cell number:','Enter start time (24h):','Enter end time (24h):','Enter depth in nucleus (microns):','Enter the channel #:','Enter age of bird on day of recording(dph):'};
prompt = {'Cell number:','Enter age of bird on day of recording(dph):'};
dlg_title = 'Input for cell statistic';
num_lines = 1;
def = {'1','70'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
% if ~(isfield(handles,'pFR'))
%     cellFile.FRsignificance=[];
% else
%     cellFile.FRsignificance=handles.pFR;
% end
cellrecnums=cell2mat(handles.cell);
%cd('C:\Users\Tim\Desktop\Tim RA Recordings');
load('C:\Users\Tim\Desktop\Tim RA Recordings\allCellsNew'); %load all Cells
cellNumber = length(allCellsNew)+1; %add the new cell to the list of cells
for i=1:length(cellrecnums)
        cellrec(i)=cellrecnums(i).filenum;
end
% buffer=0.0;


counter=1;
buffer=0.03;%don't count spikes that are within 30 ms of end of bout
for i=handles.fromIndexAudio:handles.toIndexAudio
    rec=handles.filenums(i); %filenumber we are looking at
    spikenum=find(cellrec==rec);
    if ~isempty(spikenum)
        for j=1:length(handles.audioAnnot{1,i}.segType)-1
            if ((handles.audioAnnot{1,i}.segFileStartTimes(j+1))-(handles.audioAnnot{1,i}.segFileEndTimes(j))>1&& handles.audioAnnot{1,i}.segType(j)>0)

                %an event worth looking at: find the spike file associated with
                %the filenum&& handles.audioAnnot{1,i}.segType(j)~=8 &&
                %handles.audioAnnot{1,i}.segType(j)~=9

                    firstSpikeNum=find(handles.cell{1}(spikenum).spikes>handles.audioAnnot{1,i}.segFileEndTimes(j)+buffer,1,'first');
                    postSpike(counter)=handles.cell{1}(spikenum).spikes(firstSpikeNum)-handles.audioAnnot{1,i}.segFileEndTimes(j);
                    counter=counter+1;

            end
        end


        if (handles.audioAnnot{1,i}.segFileEndTimes(end)<handles.cell{1}(spikenum).spikes(end)-1 && handles.audioAnnot{1,i}.segType(end)>0)
    %&& handles.audioAnnot{1,i}.segType(end)~=8 &&
    %handles.audioAnnot{1,i}.segType(end)~=9
                firstSpikeNum=find(handles.cell{1}(spikenum).spikes>handles.audioAnnot{1,i}.segFileEndTimes(end)+buffer,1,'first');
                postSpike(counter)=handles.cell{1}(spikenum).spikes(firstSpikeNum)-handles.audioAnnot{1,i}.segFileEndTimes(end);
                counter=counter+1;


        end
    end


end  
 %medianPS30=median(postSpike);
 allCellsNew{1,cellNumber}.cellNo=str2num(answer{1});
 allCellsNew{1,cellNumber}.dph=str2num(answer{2});
 %allCellsNew{1,cellNumber}.medianPS30 =medianPS30;

if isfield(handles, 'corr') && ~isempty(handles.corr)
    allCellsNew{1,cellNumber}.corr=handles.corr;
    allCellsNew{1,cellNumber}.motifCorrNoDrug=handles.correlation;
else
    warndlg('ai no corrida');
end

allCellsNew{1,cellNumber}.motifCorrNoDrug=handles.correlation;
allCellsNew{1,cellNumber}.AverageNoDrug=handles.motifFRND;
allCellsNew{1,cellNumber}.sparseness=handles.sparseness;
allCellsNew{1,cellNumber}.spikeI=diff(handles.spiketrainsNoDrugs{1}(1:end));
% allCellsNew{1,cellNumber}.eventFR=handles.event_firing_rate;
% allCellsNew{1,cellNumber}.eventLength =handles.event_average_length;
% allCellsNew{1,cellNumber}.event_spks =handles.event_average_spk;
 allCellsNew{1,cellNumber}.motifLength=handles.motifTemplate(end)-handles.motifTemplate(1);
for i=1:length(handles.spiketrainsNoDrugs)
    tempSpikeI=diff(handles.spiketrainsNoDrugs{i}(1:end));
    allCellsNew{1,cellNumber}.spikeI=[allCellsNew{1,cellNumber}.spikeI; tempSpikeI];
end
allCellsNew{1,cellNumber}
%cd('C:\Users\Tim\Desktop\Tim RA Recordings')
save('C:\Users\Tim\Desktop\Tim RA Recordings\allCellsNew', 'allCellsNew');
%cd( 'C:\Users\Tim\Desktop\Tim RA Recordings\Pur238 Cells for Bence')


% ISI_Callback(handles.ISI, 1, handles);
% if isfield(handles,'isiDrugCells')&& ~isempty(handles.isiDrugCells)
%     cellFile.ISIDrug{1,1}=handles.isiDrugCells(1,:);
% else
%     cellFile.ISIDrug{1,1}=[];
% end
% cellFile.ISIDrug{1,2}=handles.ifrBins;
% cellFile.ISINoDrug{1,1}=handles.isiNoDrugCells(1,:);
% cellFile.ISINoDrug{1,2}=handles.ifrBins;
% cellFile.activeDrug{1,1}=handles.activeDrug(1,:);
% cellFile.activeDrug{1,2}=handles.ifrBins;
% cellFile.activeNoDrug{1,1}=handles.activeNoDrug(1,:);
% cellFile.activeNoDrug{1,2}=handles.ifrBins;
% cellFile.frDrug=handles.frDrug;
% % cellFile.frDrug{1,2}=handles.frTimeBins;
% cellFile.PSTHcorr=handles.PSTHcorr;
% cellFile.ACF=handles.ACF;
% if isfield(handles, 'corr') && ~isempty(handles.corr)
%     cellFile.corr=handles.corr;
% else
%     warndlg('ai no corrida');
% end


% if (~isempty(handles.frNoDrug))
%     cellFile.frNoDrug=handles.frNoDrug;
% else
%     cellFile.frNoDrug=[];
% end
% cellFile.spikeI=diff(handles.spiketrainsNoDrugs{1}(1:end-1));
% cellFile.spikeI_plus=diff(handles.spiketrainsNoDrugs{1}(2:end));
% for i=1:length(handles.spiketrainsNoDrugs)
%     tempSpikeI=diff(handles.spiketrainsNoDrugs{i}(1:end-1));
%     cellFile.spikeI=[cellFile.spikeI; tempSpikeI];
%     tempSpikeI_plus=diff(handles.spiketrainsNoDrugs{i}(2:end));
%     cellFile.spikeI_plus=[cellFile.spikeI_plus; tempSpikeI_plus];
% end
% for i=1:length(cellFile.spikeI)
%     if (cellFile.spikeI(i)<0.01 && cellFile.spikeI(i)
%saves all the isi in order.

% burst=0;
% for i=1:length(cellFile.spikeI)-1
%     if (cellFile.spikeI(i)<0.0067 && cellFile.spikeI(i+1)<(0.0067))
%         burst=burst+1;
%     end
% end
% fraction=burst/length(cellFile.spikeI)
% if isempty(handles.spiketrainsDrugs)
%     cellFile.spikeIDrugs=[];
% else
%     cellFile.spikeIDrugs=diff(handles.spiketrainsDrugs{1}(1:end));
% end
% for i=1:length(handles.spiketrainsDrugs)
%     tempSpikeI=diff(handles.spiketrainsDrugs{i}(1:end));
%     cellFile.spikeIDrugs=[cellFile.spikeIDrugs; tempSpikeI];
% end
% burst=0;
% for i=1:length(cellFile.spikeIDrugs)-1
%     if (cellFile.spikeIDrugs(i)<0.0067 && cellFile.spikeIDrugs(i+1)<(0.0067))
%         burst=burst+1;
%     end
% end
% fractionDrugs=burst/length(cellFile.spikeIDrugs)
% for i=1:50 bins(i)=(0.0002)*1.15^i;
% end
% ISIHist=histc(cellFile.spikeI,bins)/length(cellFile.spikeI);
% ISIHistDrugs=histc(cellFile.spikeIDrugs,bins)/length(cellFile.spikeIDrugs);

%semilogx(bins,ISIHist);

%semilogx(bins,ISIHistDrugs,'r');
% ISIHistCum = cumsum(ISIHist);
% ISIHistDrugsCum=cumsum(ISIHistDrugs);
% figure;
% semilogx(bins,ISIHistDrugsCum,'k');
% hold on;
% semilogx(bins,ISIHistCum,'c');
% hold off;
% figure;
% loglog(bins,ISIHistDrugsCum,'k');
% hold on;
% loglog(bins,ISIHistCum,'c');
% hold off;
% figure; plot(cellFile.spikeI(1:end-1),cellFile.spikeI(2:end),'.');




    
%     hist=autocorr1(handles.spiketrainsNoDrugs{i}, res, hist_length);
%     if (i==1)
%         cellFile.stNoDrug{1,1}=hist;
%     else
%         cellFile.stNoDrug{1,1}=cellFile.stNoDrug{1,1}+hist;
%     end
% end
% cellFile.stNoDrug{1,2}=0:res:hist_length;
% for i=1:length(handles.spiketrainsDrugs)
%     hist=autocorr1(handles.spiketrainsDrugs{i}, res, hist_length);
%     if (i==1)
%         cellFile.stDrug{1,1}=hist;
%     else
%         cellFile.stDrug{1,1}=cellFile.stDrug{1,1}+hist;
%     end
% end
% 
% cellFile.stDrug{1,2}=0:res:hist_length;
%calculate autocorrelation (drug/no_drug)

% kurtosis of isi
%Isi 1 vs isi i+1
%Correlation of psth (of cells with enough repeats) with and without drugs.
%Cumulative log (isi) distribution. Logscale

% for i=1:length(handles.cell{1})
%     if (~isempty(handles.cell{1}(1,i).spont))
%         drug=find(handles.drugFiles==handles.cell{1}(1,i).filenum);
%         if isempty(drug)
%             spontSpikesNoDrug(end+1:end+length(handles.cell{1}(1,i).spont))=handles.cell{1}(1,i).spont;
%             
%         else   
%             spontSpikesDrug(end+1:end+length(handles.cell{1}(1,i).spont))=handles.cell{1}(1,i).spont;
%             
%         end
%     end
% end
% 
% cellFile.spontSpikesNoDrug=spontSpikesNoDrug;
% cellFile.spontSpikesDrug=spontSpikesDrug;
% [FileName,PathName]=uiputfile('*.mat','Save Cell File As:','Statistic_Cell');
% save([PathName FileName],'cellFile');
% showCellFile(cellFile);



% --- Executes on button press in saveSound.
function saveSound_Callback(hObject, eventdata, handles)
% hObject    handle to saveSound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%handles.fs=handles.fs;
numSyll=length(handles.sequence);
motif='';
for numSyll=1:length(handles.sequence)
    motif=[motif,num2str(handles.sequence(numSyll))]
end
%for numSyll=1:1%length(handles.sequence)
    old_dir=cd;
    [file,path]=uiputfile('',['Choose folder to save motif ',motif],'motif');
    cd (path);
    for i=1:length(handles.chosenStartSeq)
        cd(handles.dirpath);
        %open the file with the sound
        soundFileName=handles.keys{handles.chosenStartSeq(i,1)};
        [HWChannels, sound] = daq_readDatafileBence(soundFileName);
        startTime=handles.audioAnnot{1,handles.chosenStartSeq(i,1)}.segFileStartTimes(handles.chosenStartSeq(i,2))-0.01;
        endTime=handles.audioAnnot{1,handles.chosenStartSeq(i,1)}.segFileEndTimes(handles.chosenStartSeq(i,2)+(numSyll-1))+0.01;
        sylltosave=sound(round(startTime*handles.fs):round(endTime*handles.fs));
        y = resample(sylltosave,44100,handles.fs);
        cd (path);
        if (i>9)
            add='9';
        else
            add='';
        end
        wavwrite(y,44100,[num2str(handles.audioAnnot{1,handles.chosenStartSeq(i,1)}.filenum),'_motif',motif,'_',add,num2str(i)]);
        cd(handles.dirpath);

    end

%end


% --- Executes on button press in rankCorr.
function rankCorr_Callback(hObject, eventdata, handles)
% hObject    handle to rankCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


startTime=handles.motifTemplate(1)-0.06;
endTime=handles.motifTemplate(end)+0.02;
sampleRate=handles.fs;
sampRateOut=500;
GaussWidth=0.008; %in seconds
[spikeCorr,motifShifts] = crossCorrelateBence(1,handles.noMotifs, handles.spikes,sampleRate, 0, 0, GaussWidth, startTime,endTime,sampRateOut);
%[value,index]=sortMatrix(maxCorr);
%figure;imagesc(maxCorr);
[FileName,PathName]=uiputfile('*.mat','Save spike correlation As:','spikeCorr');
save([PathName FileName],'spikeCorr');

% --- Executes on button press in rankCorr.
function suppression_Callback(hObject, eventdata, handles)
% hObject    handle to suppression (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cellrecnums=cell2mat(handles.cell);

cd('C:\Users\Tim\Desktop\Tim RA Recordings');
load('allCellsNew'); %load all Cells
cellNumber = length(allCellsNew)+1;
for i=1:length(cellrecnums)
        cellrec(i)=cellrecnums(i).filenum;
end
% buffer=0.0;

for kk=1:5
    counter=1;
    buffer=(kk-1)*0.01;
    for i=handles.fromIndexAudio:handles.toIndexAudio
        rec=handles.filenums(i); %filenumber we are looking at.Type
        spikenum=find(cellrec==rec);
        if ~isempty(spikenum)
            for j=1:length(handles.audioAnnot{1,i}.segType)-1
                if ((handles.audioAnnot{1,i}.segFileStartTimes(j+1))-(handles.audioAnnot{1,i}.segFileEndTimes(j))>0.5 && handles.audioAnnot{1,i}.segType(j)>0)
                   
                    %an event worth looking at: find the spike file associated with
                    %the filenum&& handles.audioAnnot{1,i}.segType(j)~=8 &&
                    %handles.audioAnnot{1,i}.segType(j)~=9

                        firstSpikeNum=find(handles.cell{1}(spikenum).spikes>handles.audioAnnot{1,i}.segFileEndTimes(j)+buffer,1,'first');
                        postSpike(counter)=handles.cell{1}(spikenum).spikes(firstSpikeNum)-handles.audioAnnot{1,i}.segFileEndTimes(j);
                        counter=counter+1;
                    
                end
            end


            if (handles.audioAnnot{1,i}.segFileEndTimes(end)<handles.cell{1}(spikenum).spikes(end)-0.5 && handles.audioAnnot{1,i}.segType(end)>0)
        %&& handles.audioAnnot{1,i}.segType(end)~=8 &&
        %handles.audioAnnot{1,i}.segType(end)~=9
                    firstSpikeNum=find(handles.cell{1}(spikenum).spikes>handles.audioAnnot{1,i}.segFileEndTimes(end)+buffer,1,'first');
                    postSpike(counter)=handles.cell{1}(spikenum).spikes(firstSpikeNum)-handles.audioAnnot{1,i}.segFileEndTimes(end);
                    counter=counter+1;
                
                
            end
        end


    end  
    handles.medianPS(kk)=median(postSpike)
end
allCellsNew{1,cellNumber}.medianPS =handles.medianPS(1);
allCellsNew{1,cellNumber}.medianPS20 =handles.medianPS(3);
allCellsNew{1,cellNumber}.medianPS30 =handlesmedianPS(4);
% allCellsNew{1,cellNumber}.meanPS=mean(postSpike);
% allCellsNew{1,cellNumber}.stdPS=std(postSpike);
% allCellsNew{1,cellNumber}.number=counter;
allCellsNew{1,cellNumber}.number=handles.correlation;

cd('C:\Users\Tim\Desktop\Tim RA Recordings')
save('allCellsNew', 'allCellsNew');
cd('C:\Users\Tim\Desktop\Tim RA Recordings')

function cb_raster_click(hObject, evnt)
handles = guidata(hObject);
mouseMode = get(get(hObject,'Parent'), 'SelectionType');
click = get(handles.raster,'CurrentPoint');
if(strcmp(mouseMode, 'open'))
    %if double click then show the record number at click
    record=handles.audioAnnot{handles.chosenStartSeq(floor(click(1,2)))}.filenum;
    warndlg(['the record number is:',num2str(record)]);
    uiwait;
elseif(strcmp(mouseMode, 'normal') && length(handles.cell)==1)
    %zoom 
    rect = rbbox;
    endPoint = get(gca,'CurrentPoint'); 
    point1 = click(1,1:2);              % extract x and y
    point2 = endPoint(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    offset = abs(point1-point2);         % and dimensions
    l = xlim;
    if((offset(1) / l(2))< .001)%if small drag do nothing
       return;
    else
        subplot(handles.raster);
        if (~isempty(handles.eventLine1) && exist('handles.eventLine1'))
            delete(handles.eventLine1);
            delete(handles.eventLine2);
        end
        
        startTime = p1(1);
        endTime = p1(1) + offset(1);
        handles.eventLine1=line([startTime startTime],ylim, 'color',[0.2 0 0],'LineWidth',0.1);
        handles.eventLine2=line([endTime endTime],ylim, 'color',[0.2 0 0],'LineWidth',0.1);
       
        for i=1:length(handles.drugIndex)
            spikeIndex=find(handles.spikes{1,handles.drugIndex(i)}>=startTime & handles.spikes{1,handles.drugIndex(i)}<=endTime);
            drugSpikes(i)=length(handles.spikes{1,handles.drugIndex(i)}(spikeIndex));
        %calculate number of spikes in drug and no drug, average firing rate of burst in the two conditions 
        end
        for i=1:length(handles.nodrugIndex)
            spikeIndex=find(handles.spikes{1,handles.nodrugIndex(i)}>=startTime & handles.spikes{1,handles.nodrugIndex(i)}<=endTime);
            nodrugSpikes(i)=length(handles.spikes{1,handles.nodrugIndex(i)}(spikeIndex));
        %calculate number of spikes in drug and no drug, average firing rate of burst in the two conditions 
        end
        x=(0:1:max(drugSpikes));
        drugHist=histc(drugSpikes,x);
        drugHist=drugHist/sum(drugHist);
        figure;plot(drugHist,'color','r');
        hold on;
        CVDrug=std(drugSpikes)/mean(drugSpikes)
        x=(0:1:max(nodrugSpikes));
        nodrugHist=histc(nodrugSpikes,x);
        nodrugHist=nodrugHist/sum(nodrugHist);
        plot(nodrugHist);
        CVnoDrug=std(nodrugSpikes)/mean(nodrugSpikes)
        hold off;
        
        for i=1:handles.noMotifs
            spikeIndex=find(handles.spikes{1,i}>=startTime & handles.spikes{1,i}<=endTime);
            no_Spikes(i)=length(handles.spikes{1,i}(spikeIndex));
        %calculate number of spikes in drug and no drug, average firing rate of burst in the two conditions 
        end
        figure;plot(no_Spikes);
    end
        
end
guidata(hObject, handles);

% --- Executes on button press in threshold.
function threshold_Callback(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% can be found in simulateRAodyssey.m but
%more complex since it deals with both drug and nodrug recordings, and it
%has to save more information than codes that study spike trains coming
%from simulations. It is coded in an ugly way since for both drug and
%nodrug the same lines are often just copy/paste and adapted.

gaussWidth=0.003;
deriv_lim=9;
perc=8;
ratio=0.25;
rate_min=50;
diff=120;

drug_boundaries=[];
nodrug_boundaries=[];
count_previous=-1;
    if sum(handles.drugIndex==1)==1
        drug_boundaries = cat(2,drug_boundaries,1);
        count_previous=1;
    else
        nodrug_boundaries = cat(2,nodrug_boundaries,1);
        count_previous=0;
    end
    
    for i=2:handles.noMotifs
        count=sum(handles.drugIndex==i);
        if count~=count_previous & count==1
            drug_boundaries=cat(2,drug_boundaries,i);
            nodrug_boundaries=cat(2,nodrug_boundaries,i-1);
        elseif count~=count_previous & count==0
            drug_boundaries=cat(2,drug_boundaries,i-1);
            nodrug_boundaries=cat(2,nodrug_boundaries,i);
        end
        count_previous=count;
    end
    
    if count_previous==1
        drug_boundaries=cat(2,drug_boundaries,handles.noMotifs);
    else
        nodrug_boundaries=cat(2,nodrug_boundaries,handles.noMotifs);
    end

    
    
%computation of the smooth psth without FFT...

% gaussWidth=0.003;
sigma = gaussWidth / sqrt(2);
handles.bin=[];
bin_size=0.1; %bin size for psth in ms
count_index=0;
handles.psth_nodrug=[];
motif_length=(handles.motifTemplate(end)-handles.motifTemplate(1)+0.8)*1000; %length in ms
 for i=1:floor(motif_length/bin_size)
     handles.bin(i)=handles.motifTemplate(1)-0.4+(i*bin_size)/1000;
 end
 bin_size=handles.bin(2)-handles.bin(1);
 handles.psth_nodrug=zeros(1,length(handles.bin));
 [drugOut drugIn drug nodrug]=getDrugInfo(handles);
 [cND, handles.nodrugIndex] = ismember((handles.chosenStartSeq(:,1)), nodrug','rows');
 ND=find(cND==1);
 handles.nodrugIndex=ND;
 nodrugSpikes=cat(1,handles.spikes{1,handles.nodrugIndex});
 handles.psth_nodrug = hist(nodrugSpikes,handles.bin)/(length(ND)*bin_size);
handles.psth_drug=zeros(1,length(handles.bin));
 [cD, handles.drugIndex] = ismember((handles.chosenStartSeq(:,1)), drug','rows');
 D=find(cD==1);
 handles.drugIndex=D;
 drugSpikes=cat(1,handles.spikes{1,handles.drugIndex});
 handles.psth_drug = hist(drugSpikes,handles.bin)/(length(D)*bin_size);
 handles.psth_nodrug = hist(nodrugSpikes,handles.bin)/(length(ND)*bin_size);

bin_size=0.1;
for i=1:floor(motif_length/bin_size)
     handles.bin(i)=handles.motifTemplate(1)-0.4+(i*bin_size)/1000;
     if handles.bin(i)>handles.motifTemplate(1)-0.05 & count_index==0
         count_index=1;        
         index_begin=i;
     elseif  handles.bin(i)>handles.motifTemplate(end) & count_index==1
         count_index=2;
         index_end=i;
     end
end

k=0;
coeff=1;
while k*bin_size/1000 < 4*gaussWidth
    k=k+1;
    coeff=cat(2,coeff,exp(-(k*bin_size/1000)^2/(2*sigma^2)));
    coeff=cat(2,exp(-(k*bin_size/1000)^2/(2*sigma^2)),coeff);
end


for i=1:floor(length(coeff)/2) 
    coeff_start=floor(length(coeff)/2)-i+2;
    psth_start=1;
    psth_end=i+floor(length(coeff)/2);
    norm=sum(coeff(coeff_start:end));
    smooth_psth_drug(i)=sum(handles.psth_drug(psth_start:psth_end).*coeff(coeff_start:end))/norm;
    smooth_psth_nodrug(i)=sum(handles.psth_nodrug(psth_start:psth_end).*coeff(coeff_start:end))/norm;
end

for i=floor(length(coeff)/2)+1:(length(handles.bin)-floor(length(coeff)/2))

    psth_start=i-floor(length(coeff)/2);
    psth_end=i+floor(length(coeff)/2);
    smooth_psth_drug(i)=sum(handles.psth_drug(psth_start:psth_end).*coeff)/norm;
    smooth_psth_nodrug(i)=sum(handles.psth_nodrug(psth_start:psth_end).*coeff)/norm;
end

for i=((length(handles.bin)-floor(length(coeff)/2))+1):length(handles.bin)
    coeff_end=length(coeff)-i+length(handles.bin)-floor(length(coeff)/2);
    psth_start=i-floor(length(coeff)/2);
    psth_end=length(handles.bin);
    smooth_psth_drug(i)=sum(handles.psth_drug(psth_start:psth_end).*coeff(1:coeff_end))/norm;
    smooth_psth_nodrug(i)=sum(handles.psth_nodrug(psth_start:psth_end).*coeff(1:coeff_end))/norm;
end


%compute the derivative of the smooth psth
% deriv_lim=9;
deriv_psth_drug = zeros(floor(motif_length/bin_size),1);
for i=2:floor(motif_length/bin_size)
    deriv_psth_drug(i) = (smooth_psth_drug(i)-smooth_psth_drug(i-1))/bin_size;
end

deriv_psth_nodrug = zeros(floor(motif_length/bin_size),1);
for i=2:floor(motif_length/bin_size)
    deriv_psth_nodrug(i) = (smooth_psth_nodrug(i)-smooth_psth_nodrug(i-1))/bin_size;
end
max_deriv=max(abs(deriv_psth_nodrug(index_begin:index_end)));
if max_deriv<deriv_lim
    seq=num2str(handles.Seq);
    nb_drug_renditions=num2str(length(handles.drugIndex));
    spk=handles.spikes;
    motifTemplate_start=handles.motifTemplate(1);
    motifTemplate_end=handles.motifTemplate(end);
    file_name_nodrug=strcat('nodrug_',handles.name_cell,'_syll',seq,'_deriv');
    save(file_name_nodrug,'spk','motifTemplate_start','motifTemplate_end')


elseif length(handles.nodrugIndex)<11
    seq=num2str(handles.Seq);
    nb_drug_renditions=num2str(length(handles.drugIndex));
    spk=handles.spikes;
    motifTemplate_start=handles.motifTemplate(1);
    motifTemplate_end=handles.motifTemplate(end);
    file_name_nodrug=strcat('nodrug_',handles.name_cell,'_syll',seq,'_nbrend');
    save(file_name_nodrug,'spk','motifTemplate_start','motifTemplate_end')
    
    
else


%spot the maximum and the minimum of the psth.

%Variables that start with
%the word time code for a time in the sequence (final variable). Variables
%that start with the word index code for the index (relative to the vector
%coding for smooth psth) of those time (intermediary variable)

 max_drug=[];
 min_drug=[0,0];
 max_nodrug=[];
 min_nodrug=[0,0];
 index_max_drug=[];
 index_min_drug=[0,0];
 index_max_nodrug=[];
 index_min_nodrug=[0,0];
 time_max_drug=[];
 time_min_drug=[handles.motifTemplate(1)-0.1,handles.motifTemplate(end)+0.1];
 time_max_nodrug=[];
 time_min_nodrug=[handles.motifTemplate(1)-0.1,handles.motifTemplate(end)+0.1];
 nodrug_time_events_l=[];
 nodrug_time_events_r=[];
 drug_time_events_l=[];
 drug_time_events_r=[];
 nodrug_index_events_l=[];
 nodrug_index_events_r=[];
 drug_index_events_l=[];
 drug_index_events_r=[];
 count_drug=0;
 count_drug_max=0;
 count_nodrug=0;
 count_nodrug_max=0;
 t_prev_max_drug=-100;
t_prev_min_drug=-10;
t_prev_max_nodrug=-100;
t_prev_min_nodrug=-10;

for i=2:length(handles.bin)-1
    if handles.bin(i)>handles.motifTemplate(1)-0.1 & handles.bin(i)<handles.motifTemplate(end)+0.1
    if smooth_psth_drug(i+1)>smooth_psth_drug(i) & smooth_psth_drug(i-1)>smooth_psth_drug(i)
        if count_drug_max~=0
        count_drug=count_drug-1;
        end
       min_drug=cat(2,min_drug,smooth_psth_drug(i));
       time_min_drug=cat(2,time_min_drug,handles.bin(i));
       index_min_drug=cat(2,index_min_drug,i);
       t_prev_min_drug=handles.bin(i);

    elseif smooth_psth_drug(i+1)<smooth_psth_drug(i) & smooth_psth_drug(i-1)<smooth_psth_drug(i) %& smooth_psth_drug(i)>0.5
       count_drug=count_drug+1;
       count_drug_max=count_drug_max+1;
       max_drug=cat(2,max_drug,smooth_psth_drug(i));
       time_max_drug=cat(2,time_max_drug,handles.bin(i));
       index_max_drug=cat(2,index_max_drug,i);
       if count_drug==2
           if length(index_max_drug)>1
           [min_ad index_min_ad]=min(smooth_psth_drug(index_max_drug(end-1):index_max_drug(end)));
           min_drug=cat(2,min_drug,smooth_psth_drug(index_min_ad+index_max_drug(end-1)-1));
           time_min_drug=cat(2,time_min_drug,handles.bin(index_min_ad+index_max_drug(end-1)-1));
           index_min_drug=cat(2,index_min_drug,index_min_ad+index_max_drug(end-1)-1);
           count_drug=0;
           end
       end
       t_prev_max_drug=handles.bin(i);
       elseif smooth_psth_drug(i-1)==smooth_psth_drug(i) & smooth_psth_drug(i)>5
        k=smooth_psth_drug(i);
        l=i;
        while k==smooth_psth_drug(i)
            l=l+1;
            k=smooth_psth_drug(l);
        end
        m=smooth_psth_drug(i-1);
        n=i-1;
        while m==smooth_psth_drug(i-1)
            n=n-1;
            m=smooth_psth_drug(n);
        end
        if (m-smooth_psth_drug(i-1))<0 & (k-smooth_psth_drug(i-1))<0
        
         count_drug=count_drug+1;
       count_drug_max=count_drug_max+1;
       max_drug=cat(2,max_drug,smooth_psth_drug(i));
       time_max_drug=cat(2,time_max_drug,handles.bin(i));
       index_max_drug=cat(2,index_max_drug,i);
       if count_drug==2
  
           [min_ad index_min_ad]=min(smooth_psth_drug(index_max_drug(end-1):index_max_drug(end)));
           min_drug=cat(2,min_drug,smooth_psth_drug(index_min_ad+index_max_drug(end-1)-1));
           time_min_drug=cat(2,time_min_drug,handles.bin(index_min_ad+index_max_drug(end-1)-1));
           index_min_drug=cat(2,index_min_drug,index_min_ad+index_max_drug(end-1)-1);
       end
        t_prev_max_drug=handles.bin(i);
        end
    end
    
    %loop to avoid to spot min and max that are too close from each other (cf low rate events)
    if abs(t_prev_max_drug-t_prev_min_drug)<gaussWidth/2
        time_min_drug=time_min_drug(1:end-1);
        index_min_drug=index_min_drug(1:end-1);
        min_drug=min_drug(1:end-1);
        time_max_drug=time_max_drug(1:end-1);
        index_max_drug=index_max_drug(1:end-1);
        max_drug=max_drug(1:end-1);
        t_prev_max_drug=-100;
        t_prev_min_drug=-10;
    end

    if smooth_psth_nodrug(i+1)>smooth_psth_nodrug(i) & smooth_psth_nodrug(i-1)>smooth_psth_nodrug(i)
       
        if count_nodrug_max~=0
        count_nodrug=count_nodrug-1;
        end
        min_nodrug=cat(2,min_nodrug,smooth_psth_nodrug(i));
       time_min_nodrug=cat(2,time_min_nodrug,handles.bin(i));
       index_min_nodrug=cat(2,index_min_nodrug,i);
       t_prev_min_nodrug=handles.bin(i);

    elseif smooth_psth_nodrug(i+1)<smooth_psth_nodrug(i) & smooth_psth_nodrug(i-1)<smooth_psth_nodrug(i)%& smooth_psth_drug(i)>0.5
         
        count_nodrug=count_nodrug+1;
       count_nodrug_max=count_nodrug_max+1;
       max_nodrug=cat(2,max_nodrug,smooth_psth_nodrug(i));
       time_max_nodrug=cat(2,time_max_nodrug,handles.bin(i));
       index_max_nodrug=cat(2,index_max_nodrug,i);
       t_prev_max_nodrug=handles.bin(i);
       %the folowing loop ad a minimum between two maximums when the code
       %is unable to spot this minimum
       if count_nodrug==2
           if length(index_max_nodrug)>1
           [min_ad index_min_ad]=min(smooth_psth_nodrug(index_max_nodrug(end-1):index_max_nodrug(end)));
           min_nodrug=cat(2,min_nodrug,smooth_psth_nodrug(index_min_ad+index_max_nodrug(end-1)-1));
           time_min_nodrug=cat(2,time_min_nodrug,handles.bin(index_min_ad+index_max_nodrug(end-1)-1));
           index_min_nodrug=cat(2,index_min_nodrug,index_min_ad+index_max_nodrug(end-1)-1);
           end
       end
    elseif smooth_psth_nodrug(i-1)==smooth_psth_nodrug(i) & smooth_psth_nodrug(i)>5
        k=smooth_psth_nodrug(i);
        l=i;
        while k==smooth_psth_nodrug(i)
            l=l+1;
            k=smooth_psth_nodrug(l);
        end
        m=smooth_psth_nodrug(i-1);
        n=i-1;
        while m==smooth_psth_nodrug(i-1)
            n=n-1;
            m=smooth_psth_nodrug(n);
        end
        if (m-smooth_psth_nodrug(i-1))<0 & (k-smooth_psth_nodrug(i-1))<0
        
         count_nodrug=count_nodrug+1;
       count_nodrug_max=count_nodrug_max+1;
       max_nodrug=cat(2,max_nodrug,smooth_psth_nodrug(i));
       time_max_nodrug=cat(2,time_max_nodrug,handles.bin(i));
       index_max_nodrug=cat(2,index_max_nodrug,i);
       if count_nodrug==2
           [min_ad index_min_ad]=min(smooth_psth_nodrug(index_max_nodrug(end-1):index_max_nodrug(end)));
           min_nodrug=cat(2,min_nodrug,smooth_psth_nodrug(index_min_ad+index_max_nodrug(end-1)-1));
           time_min_nodrug=cat(2,time_min_nodrug,handles.bin(index_min_ad+index_max_nodrug(end-1)-1));
           index_min_nodrug=cat(2,index_min_nodrug,index_min_ad+index_max_nodrug(end-1)-1);
       end
        t_prev_max_nodrug=handles.bin(i);
        end
    end
    if abs(t_prev_max_nodrug-t_prev_min_nodrug)<gaussWidth/2
        time_min_nodrug=time_min_nodrug(1:end-1);
        index_min_nodrug=index_min_nodrug(1:end-1);
        min_nodrug=min_nodrug(1:end-1);
        time_max_nodrug=time_max_nodrug(1:end-1);
        index_max_nodrug=index_max_nodrug(1:end-1);
        max_nodrug=max_nodrug(1:end-1);
        t_prev_max_nodrug=-100;
        t_prev_min_nodrug=-10;
    end
    end
end



%compute the time where the smooth psth cross the percentage of the maximum 
num_nodrug=[];
temp_max=max_nodrug;
temp_index=index_max_nodrug;
nodrug_time_events=[];
while length(temp_max)>0
    
     num_nodrug=num_nodrug+1;  %give the number of the event
     index_real=temp_index(1);
     maxi=temp_max(1);
     temp_max=temp_max(2:end);
     temp_index=temp_index(2:end);
    
    
%     perc=10; %percentage of the studied maximum coresponding to the time
    %that can be the boundaries of the events
    k=index_real;
    p=smooth_psth_nodrug(index_real);
   
    while p > perc*maxi/100 & handles.bin(k+1)<handles.motifTemplate(end)+0.1
        k=k+1;
        p = smooth_psth_nodrug(k);
    end
        q=smooth_psth_nodrug(index_real);
        l=index_real;
    while q > perc*maxi/100 & handles.bin(l-1)>handles.motifTemplate(1)-0.1
        l=l-1;
        q = smooth_psth_nodrug(l);
    end
    %spot the left and right minimum of this maximum
    %handles.bin(index_real)= time of the studied maximum
    min_ind = time_min_nodrug-handles.bin(index_real);
    for i=1:length(min_ind)
        if min_ind(i)>=0
            min_ind_right(i)=min_ind(i);
            min_ind_left(i)=-100;
        else
            min_ind_right(i)=100;
            min_ind_left(i)=min_ind(i);
        end
    end
        [m_right i_right]=min(min_ind_right);
        [m_left i_left]=max(min_ind_left);
        t_min_right=time_min_nodrug(i_right);
        t_min_left=time_min_nodrug(i_left);
    
    if handles.bin(l)> t_min_left
        nodrug_time_events_l=cat(2,nodrug_time_events_l,handles.bin(l));
        nodrug_index_events_l=cat(2,nodrug_index_events_l,l);
    else
        nodrug_time_events_l=cat(2,nodrug_time_events_l,t_min_left);
        nodrug_index_events_l=cat(2,nodrug_index_events_l,index_min_nodrug(i_left));
        
    end
    if handles.bin(k)<t_min_right
        nodrug_time_events_r=cat(2,nodrug_time_events_r,handles.bin(k));
        nodrug_index_events_r=cat(2,nodrug_index_events_r,k);
    else
        nodrug_time_events_r=cat(2,nodrug_time_events_r,t_min_right);
        nodrug_index_events_r=cat(2,nodrug_index_events_r,index_min_nodrug(i_right));
    end
    
end



num_drug=0;
temp_max=max_drug;
temp_index=index_max_drug;
drug_time_events=[];
while length(temp_max)>0
    
     num_drug=num_drug+1;  %give the number of the event
     index_real=temp_index(1);
     maxi=temp_max(1);
     temp_max=temp_max(2:end);
     temp_index=temp_index(2:end);

%     perc=10;
    k=index_real;
    p=smooth_psth_drug(index_real);
    while p > perc*maxi/100 & handles.bin(k+1)<handles.motifTemplate(end)+0.1
        k=k+1;
        p = smooth_psth_drug(k);
    end
        q=smooth_psth_drug(index_real);
        l=index_real;
    while q > perc*maxi/100 & handles.bin(l-1)>handles.motifTemplate(1)-0.1
        l=l-1;
        q = smooth_psth_drug(l);
    end
    
    %spot the left and right minimum of this maximum
    %handles.bin(index_real)=time of the studied maximum
    min_ind = time_min_drug-handles.bin(index_real);
    min_ind_right=[];
    min_ind_left=[];
    for i=1:length(min_ind)
        if min_ind(i)>=0
            min_ind_right(i)=min_ind(i);
            min_ind_left(i)=-100;
        else
            min_ind_right(i)=100;
            min_ind_left(i)=min_ind(i);
        end
    end

        [m_right i_right]=min(min_ind_right);
        [m_left i_left]=max(min_ind_left);
        t_min_right=time_min_drug(i_right);
        t_min_left=time_min_drug(i_left);

    if handles.bin(l)> t_min_left
        drug_time_events_l=cat(2,drug_time_events_l,handles.bin(l));
        drug_index_events_l=cat(2,drug_index_events_l,l);
    else
        drug_time_events_l=cat(2,drug_time_events_l,t_min_left);
        drug_index_events_l=cat(2,drug_index_events_l,index_min_drug(i_left));
    end
    if handles.bin(k)<t_min_right
        drug_time_events_r=cat(2,drug_time_events_r,handles.bin(k));
        drug_index_events_r=cat(2,drug_index_events_r,k);
    else
        drug_time_events_r=cat(2,drug_time_events_r,t_min_right);
        drug_index_events_r=cat(2,drug_index_events_r,index_min_drug(i_right));
    end
    
end



%loops to reject events with not enough spikes
% ratio=0.25; %if nb_of_spikes/nb_of_renditions<ratio then reject event
keep_index=[];
for k=1:length(nodrug_time_events_r)
    t_inf=nodrug_time_events_l(k);
    t_sup=nodrug_time_events_r(k);
  for j=1:length(handles.cell)
        for i=1:handles.noMotifs
            spikeIndex=find(handles.spikes{j,i}>=t_inf & handles.spikes{j,i}<=t_sup);
            test = handles.nodrugIndex == i;
            if sum(test)==0
            handles.nodrug_spikes{j,i}=[];
            else
            handles.nodrug_spikes{j,i}=handles.spikes{j,i}(spikeIndex);
            end
        end
  end
  nb_of_spikes=0;
  for i=1:handles.noMotifs
      nb_of_spikes=nb_of_spikes+length(handles.nodrug_spikes{i});
  end
  if nb_of_spikes/handles.noMotifs>ratio
      keep_index=cat(2,keep_index,k);
  end
end
nodrug_time_events_l=nodrug_time_events_l(keep_index);
nodrug_time_events_r=nodrug_time_events_r(keep_index);
nodrug_index_events_l=nodrug_index_events_l(keep_index);
nodrug_index_events_r=nodrug_index_events_r(keep_index);




% %reject events out of bounds


keep_index=[];
for k=1:length(drug_time_events_r)
    t_inf=drug_time_events_l(k);
    t_sup=drug_time_events_r(k);
  for j=1:length(handles.cell)
        for i=1:handles.noMotifs
            spikeIndex=find(handles.spikes{j,i}>=t_inf & handles.spikes{j,i}<=t_sup);
            test = handles.drugIndex == i;
            if sum(test)==0
            handles.drug_spikes{j,i}=[];
            else
            handles.drug_spikes{j,i}=handles.spikes{j,i}(spikeIndex);
            end
        end
  end
  nb_of_spikes=0;
  for i=1:handles.noMotifs
      nb_of_spikes=nb_of_spikes+length(handles.drug_spikes{i});
  end
  
  if nb_of_spikes/handles.noMotifs>ratio
      keep_index=cat(2,keep_index,k);
  end
end
drug_time_events_l=drug_time_events_l(keep_index);
drug_time_events_r=drug_time_events_r(keep_index);
drug_index_events_l=drug_index_events_l(keep_index);
drug_index_events_r=drug_index_events_r(keep_index);

% %reject events out of bounds
% keep_index=[];
% for k=1:length(drug_time_events_r)
%     if drug_time_events_l(k)>handles.motifTemplate(1)-0.05 & ...
%             drug_time_events_r(k)<handles.motifTemplate(end)
%         keep_index=cat(2,keep_index,k);
%     end
% end
% drug_time_events_l=drug_time_events_l(keep_index);
% drug_time_events_r=drug_time_events_r(keep_index);
% drug_index_events_l=drug_index_events_l(keep_index);
% drug_index_events_r=drug_index_events_r(keep_index);



%reject events whose max psth < 30Hz
% rate_min=30;
keep_index=[];
for i=1:length(nodrug_index_events_r)
    if max(smooth_psth_nodrug(nodrug_index_events_l(i):nodrug_index_events_r(i)))>rate_min
        keep_index=cat(2,keep_index,i);
    end
end
nodrug_time_events_l=nodrug_time_events_l(keep_index);
nodrug_time_events_r=nodrug_time_events_r(keep_index);
nodrug_index_events_l=nodrug_index_events_l(keep_index);
nodrug_index_events_r=nodrug_index_events_r(keep_index);

keep_index=[];
for i=1:length(drug_index_events_r)
    max(smooth_psth_drug(drug_index_events_l(i):drug_index_events_r(i)));
    smooth_psth_drug(drug_index_events_l(i));
    if max(smooth_psth_drug(drug_index_events_l(i):drug_index_events_r(i)))>rate_min
        keep_index=cat(2,keep_index,i);
    end
end
drug_time_events_l=drug_time_events_l(keep_index);
drug_time_events_r=drug_time_events_r(keep_index);
drug_index_events_l=drug_index_events_l(keep_index);
drug_index_events_r=drug_index_events_r(keep_index);




        
%merging events-------------------------
% diff=20;
nodrug_time_events_r_merged=[];
nodrug_index_events_r_merged=[];

if isempty(nodrug_time_events_l)==0
nodrug_time_events_l_merged=nodrug_time_events_l(1);
nodrug_index_events_l_merged=nodrug_index_events_l(1);
for i=2:length(nodrug_time_events_l)
    maxi1=max(smooth_psth_nodrug(nodrug_index_events_l(i-1):nodrug_index_events_r(i-1)));
    maxi2=max(smooth_psth_nodrug(nodrug_index_events_l(i):nodrug_index_events_r(i)));
        %-> try to merge events if they have a common boundary
    if nodrug_time_events_l(i)==nodrug_time_events_r(i-1)
        %-> condition for merging (see rules definition)
        if min(maxi1,maxi2)-smooth_psth_nodrug(nodrug_index_events_l(i))<diff
        else
            nodrug_time_events_l_merged=cat(2,nodrug_time_events_l_merged,nodrug_time_events_l(i));
            nodrug_time_events_r_merged=cat(2,nodrug_time_events_r_merged,nodrug_time_events_r(i-1));
            nodrug_index_events_l_merged=cat(2,nodrug_index_events_l_merged,nodrug_index_events_l(i));
            nodrug_index_events_r_merged=cat(2,nodrug_index_events_r_merged,nodrug_index_events_r(i-1));
        end
    else
        nodrug_time_events_l_merged=cat(2,nodrug_time_events_l_merged,nodrug_time_events_l(i));
        nodrug_time_events_r_merged=cat(2,nodrug_time_events_r_merged,nodrug_time_events_r(i-1));
        nodrug_index_events_l_merged=cat(2,nodrug_index_events_l_merged,nodrug_index_events_l(i));
        nodrug_index_events_r_merged=cat(2,nodrug_index_events_r_merged,nodrug_index_events_r(i-1));
        
    end
end
%-> in any case, add the last right boundary spotted previously
nodrug_time_events_r_merged=cat(2,nodrug_time_events_r_merged,nodrug_time_events_r(end));
nodrug_index_events_r_merged=cat(2,nodrug_index_events_r_merged,nodrug_index_events_r(end));
else
    nodrug_time_events_r_merged=[];
    nodrug_time_events_l_merged=[];
    nodrug_index_events_r_merged=[];
    nodrug_index_events_l_merged=[];
end
%-> same for drug
if isempty(drug_time_events_l)==0
drug_time_events_r_merged=[];
drug_time_events_l_merged=drug_time_events_l(1);
drug_index_events_r_merged=[];
drug_index_events_l_merged=drug_index_events_l(1);
for i=2:length(drug_time_events_l)
    maxi1=max(smooth_psth_drug(drug_index_events_l(i-1):drug_index_events_r(i-1)));
    maxi2=max(smooth_psth_drug(drug_index_events_l(i):drug_index_events_r(i)));
    if drug_time_events_l(i)==drug_time_events_r(i-1)
        if min(maxi1,maxi2)-smooth_psth_drug(drug_index_events_l(i))<diff
        else
            drug_time_events_l_merged=cat(2,drug_time_events_l_merged,drug_time_events_l(i));
            drug_time_events_r_merged=cat(2,drug_time_events_r_merged,drug_time_events_r(i-1));
            drug_index_events_l_merged=cat(2,drug_index_events_l_merged,drug_index_events_l(i));
            drug_index_events_r_merged=cat(2,drug_index_events_r_merged,drug_index_events_r(i-1));
        end
    else
        drug_time_events_l_merged=cat(2,drug_time_events_l_merged,drug_time_events_l(i));
        drug_time_events_r_merged=cat(2,drug_time_events_r_merged,drug_time_events_r(i-1));
        drug_index_events_l_merged=cat(2,drug_index_events_l_merged,drug_index_events_l(i));
        drug_index_events_r_merged=cat(2,drug_index_events_r_merged,drug_index_events_r(i-1));
    end
end
drug_time_events_r_merged=cat(2,drug_time_events_r_merged,drug_time_events_r(end));
drug_index_events_r_merged=cat(2,drug_index_events_r_merged,drug_index_events_r(end));
else
    drug_time_events_r_merged=[];
    drug_time_events_l_merged=[];
    drug_index_events_r_merged=[];
    drug_index_events_l_merged=[];
end



%-> redefine the boundaries of the events, given that we now have a new
%time for the maximum (last part of rule 2 reapply to those new minimum)
for i=1:length(nodrug_index_events_l_merged)
    count_mod=0;
    [maxim index_bis]=max( smooth_psth_nodrug(...
        nodrug_index_events_l_merged(i):nodrug_index_events_r_merged(i)) );
    if maxim>diff
        %-> spot new times associated to the new maximum
    k=maxim;
    l=nodrug_index_events_l_merged(i)+index_bis-1;
    while k>=perc/100*maxim & l<length(handles.bin)
        l=l+1;
        k=smooth_psth_nodrug(l);
    end
    m=maxim;
    n=nodrug_index_events_l_merged(i)+index_bis-1;
    while m>=perc/100*maxim & n>1
        n=n-1;
        handles.bin(n);
        m=smooth_psth_nodrug(n);
    end
    if l<nodrug_index_events_r_merged(i)
        count_mod=count_mod+1;
   ad_index_right=l;
   ad_time_right=handles.bin(l);
    end
     if n>nodrug_index_events_l_merged(i)
       count_mod=count_mod+0.5;
   ad_index_left=n;
   ad_time_left=handles.bin(n);
     end
     
      %-> decide which of those times have to be defined as the new events
      %boundaries according to rule 2
    if count_mod==1.5
        if nodrug_index_events_l_merged(i+1)==nodrug_index_events_r_merged(i)...
                & nodrug_index_events_l_merged(i)==nodrug_index_events_r_merged(i-1)
     nodrug_index_events_l_merged(i+1)=ad_index_right;
       nodrug_time_events_l_merged(i+1)=ad_time_right;
       nodrug_index_events_r_merged(i)=ad_index_right;
       nodrug_time_events_r_merged(i)=ad_time_right;
       nodrug_index_events_l_merged(i)=ad_index_left;
       nodrug_time_events_l_merged(i)=ad_time_left;
       nodrug_index_events_r_merged(i-1)=ad_index_left;
       nodrug_time_events_r_merged(i-1)=ad_time_left;
        elseif nodrug_index_events_l_merged(i+1)==nodrug_index_events_r_merged(i)
            nodrug_index_events_l_merged(i+1)=ad_index_right;
       nodrug_time_events_l_merged(i+1)=ad_time_right;
       nodrug_index_events_r_merged(i)=ad_index_right;
       nodrug_time_events_r_merged(i)=ad_time_right;
       nodrug_index_events_l_merged(i)=ad_index_left;
       nodrug_time_events_l_merged(i)=ad_time_left;
        elseif nodrug_index_events_l_merged(i)==nodrug_index_events_r_merged(i-1)
            nodrug_index_events_r_merged(i)=ad_index_right;
       nodrug_time_events_r_merged(i)=ad_time_right;
       nodrug_index_events_l_merged(i)=ad_index_left;
       nodrug_time_events_l_merged(i)=ad_time_left;
       nodrug_index_events_r_merged(i-1)=ad_index_left;
       nodrug_time_events_r_merged(i-1)=ad_time_left;
        else
            nodrug_index_events_r_merged(i)=ad_index_right;
       nodrug_time_events_r_merged(i)=ad_time_right;
       nodrug_index_events_l_merged(i)=ad_index_left;
       nodrug_time_events_l_merged(i)=ad_time_left;
        end
    elseif count_mod==1
        if i==length(nodrug_index_events_l_merged)+1
            nodrug_index_events_r_merged(i)=ad_index_right;
       nodrug_time_events_r_merged(i)=ad_time_right;
        else
        if (i<length(nodrug_index_events_l_merged) && nodrug_index_events_l_merged(i+1)==nodrug_index_events_r_merged(i))
       nodrug_index_events_l_merged(i+1)=ad_index_right;
       nodrug_time_events_l_merged(i+1)=ad_time_right;
       nodrug_index_events_r_merged(i)=ad_index_right;
       nodrug_time_events_r_merged(i)=ad_time_right;
        else
       nodrug_index_events_r_merged(i)=ad_index_right;
       nodrug_time_events_r_merged(i)=ad_time_right;
        end
        end
   elseif count_mod==0.5
       if i-1==0
           nodrug_index_events_l_merged(i)=ad_index_left;
           nodrug_time_events_l_merged(i)=ad_time_left;
       else
       if nodrug_index_events_l_merged(i)==nodrug_index_events_r_merged(i-1)       
           nodrug_index_events_l_merged(i)=ad_index_left;
           nodrug_time_events_l_merged(i)=ad_time_left;
           nodrug_index_events_r_merged(i-1)=ad_index_left;
           nodrug_time_events_r_merged(i-1)=ad_time_left;
       else
       nodrug_index_events_l_merged(i)=ad_index_left;
       nodrug_time_events_l_merged(i)=ad_time_left;
       end
       end
    end
    end
end

%reject events out of bounds
keep_index=[];
for k=1:length(nodrug_time_events_r_merged)
    if nodrug_time_events_l_merged(k)>handles.motifTemplate(1)-0.05 & ...
            nodrug_time_events_r_merged(k)<handles.motifTemplate(end)
        keep_index=cat(2,keep_index,k);
    end
end
nodrug_time_events_l_merged=nodrug_time_events_l_merged(keep_index);
nodrug_time_events_r_merged=nodrug_time_events_r_merged(keep_index);
nodrug_index_events_l_merged=nodrug_index_events_l_merged(keep_index);
nodrug_index_events_r_merged=nodrug_index_events_r_merged(keep_index);

%reject event with low psth derivative
keep_index=[];
for k=1:length(nodrug_time_events_r_merged)
    max_deriv=max(abs(deriv_psth_nodrug(nodrug_index_events_l_merged(k):nodrug_index_events_r_merged(k))));
if max_deriv>=deriv_lim
    keep_index=cat(2,keep_index,k);
end
end
nodrug_time_events_l_merged=nodrug_time_events_l_merged(keep_index);
nodrug_time_events_r_merged=nodrug_time_events_r_merged(keep_index);
nodrug_index_events_l_merged=nodrug_index_events_l_merged(keep_index);
nodrug_index_events_r_merged=nodrug_index_events_r_merged(keep_index);
 




for i=1:length(drug_index_events_l_merged)
    count_mod=0;
%     drug_time_events_l_merged
%     drug_time_events_r_merged
    [maxim index_bis]=max( smooth_psth_drug(...
        drug_index_events_l_merged(i):drug_index_events_r_merged(i)) );
    if maxim>diff
    k=maxim;
    l=drug_index_events_l_merged(i)+index_bis-1;
    while k>=perc/100*maxim & l<length(handles.bin)
        l=l+1;
        k=smooth_psth_drug(l);
    end
    m=maxim;
    n=drug_index_events_l_merged(i)+index_bis-1;
    while m>=perc/100*maxim & n>1
        n=n-1;
        handles.bin(n);
        m=smooth_psth_drug(n);
    end
    if l<drug_index_events_r_merged(i)
        count_mod=count_mod+1;
   ad_index_right=l;
   ad_time_right=handles.bin(l);
    end
     if n>drug_index_events_l_merged(i)
       count_mod=count_mod+0.5;
   ad_index_left=n;
   ad_time_left=handles.bin(n);
     end
     
    if count_mod==1.5
         drug_index_events_l_merged=cat(2,drug_index_events_l_merged(1:i),...
           ad_index_left,ad_index_right,drug_index_events_l_merged(i+1:end));
       drug_time_events_l_merged=cat(2,drug_time_events_l_merged(1:i),...
           ad_time_left,ad_time_right,drug_time_events_l_merged(i+1:end));
       drug_index_events_r_merged=cat(2,drug_index_events_r_merged(1:i-1),...
           ad_index_left,ad_index_right,drug_index_events_r_merged(i:end));
       drug_time_events_r_merged=cat(2,drug_time_events_r_merged(1:i-1),...
           ad_time_left,ad_time_right,drug_time_events_r_merged(i:end)) ;
    elseif count_mod==1
       drug_index_events_l_merged=cat(2,drug_index_events_l_merged(1:i),...
           ad_index_right,drug_index_events_l_merged(i+1:end));
       drug_time_events_l_merged=cat(2,drug_time_events_l_merged(1:i),...
           ad_time_right,drug_time_events_l_merged(i+1:end));
       drug_index_events_r_merged=cat(2,drug_index_events_r_merged(1:i-1),...
           ad_index_right,drug_index_events_r_merged(i:end));
       drug_time_events_r_merged=cat(2,drug_time_events_r_merged(1:i-1),...
           ad_time_right,drug_time_events_r_merged(i:end));
   elseif count_mod==0.5
       drug_index_events_l_merged=cat(2,drug_index_events_l_merged(1:i),...
           ad_index_left,drug_index_events_l_merged(i+1:end));
       drug_time_events_l_merged=cat(2,drug_time_events_l_merged(1:i),...
           ad_time_left,drug_time_events_l_merged(i+1:end));
       drug_index_events_r_merged=cat(2,drug_index_events_r_merged(1:i-1),...
           ad_index_left,drug_index_events_r_merged(i:end));
       drug_time_events_r_merged=cat(2,drug_time_events_r_merged(1:i-1),...
           ad_time_left,drug_time_events_r_merged(i:end));
    end
    end
end

%reject events drug out of bounds
keep_index=[];
for k=1:length(drug_time_events_r_merged)
    if drug_time_events_l_merged(k)>handles.motifTemplate(1)-0.05 & ...
            drug_time_events_r_merged(k)<handles.motifTemplate(end)
        keep_index=cat(2,keep_index,k);
    end
end
drug_time_events_l_merged=drug_time_events_l_merged(keep_index);
drug_time_events_r_merged=drug_time_events_r_merged(keep_index);
drug_index_events_l_merged=drug_index_events_l_merged(keep_index);
drug_index_events_r_merged=drug_index_events_r_merged(keep_index);



%reject event with low psth derivative
keep_index=[];
for k=1:length(drug_time_events_r_merged)
    max_deriv=max(abs(deriv_psth_drug(drug_index_events_l_merged(k):drug_index_events_r_merged(k))));
if max_deriv>=deriv_lim
    keep_index=cat(2,keep_index,k);
end
end
drug_time_events_l_merged=drug_time_events_l_merged(keep_index);
drug_time_events_r_merged=drug_time_events_r_merged(keep_index);
drug_index_events_l_merged=drug_index_events_l_merged(keep_index);
drug_index_events_r_merged=drug_index_events_r_merged(keep_index);





%save the events

handles.nodrug_spikes=[];
handles.drug_spikes=[];
for k=1:length(drug_index_events_r_merged)
    t_inf=drug_time_events_l_merged(k);
    t_sup=drug_index_events_r_merged(k);
    for j=1:length(handles.cell)
        for i=1:handles.noMotifs
            spikeIndex=find(handles.spikes{j,i}>=t_inf & handles.spikes{j,i}<=t_sup);
            test = handles.drugIndex == i;
            if sum(test)==0
            handles.drug_spikes{j,i}=[-100];
            handles.nodrug_spikes{j,i}=handles.spikes{j,i}(spikeIndex);
            else
            handles.drug_spikes{j,i}=handles.spikes{j,i}(spikeIndex);
            handles.nodrug_spikes{j,i}=[-100];
            end
        end
    end


[maxi index_max]=max(smooth_psth_drug(drug_index_events_l_merged(k):drug_index_events_r_merged(k)));
time_max=handles.bin(drug_index_events_l_merged(k)+index_max-1);
m_left=smooth_psth_drug(drug_index_events_l_merged(k));
m_right=smooth_psth_drug(drug_index_events_r_merged(k));
motifTemplate_start=handles.motifTemplate(1);
motifTemplate_end=handles.motifTemplate(end);
t_min_left=drug_time_events_l_merged(k);
t_min_right=drug_time_events_r_merged(k);
num=k;
num_str=num2str(num);
zero_str=num2str(0);
spk=handles.drug_spikes;
spk_bis=handles.nodrug_spikes;
%gives a status to the event
nb_of_max = t_min_left<time_max_drug  & time_max_drug<t_min_right;
if sum(nb_of_max)==1
    status=1;
elseif sum(nb_of_max)==2 & maxi>20
    status=2; %will identify the double bump events
else
    status=3;%should correspond to low rate events
end
    

if num<10
    num_str=strcat(zero_str,num_str);
end
seq=num2str(handles.Seq);
nb_drug_renditions=num2str(length(handles.drugIndex));
file_name_drug=strcat('drug_',handles.name_cell,'_syll',seq,'_event',num_str);
drugIndex=handles.drugIndex;
% if num~=(length(time_min_drug)-1)
%     plant
% else
save(file_name_drug,'spk','spk_bis','nb_drug_renditions','t_min_left','t_min_right'...
    ,'time_max','m_left','m_right','maxi','motifTemplate_start','motifTemplate_end',...
    'status','drugIndex')
%end
end

handles.nodrug_spikes=[];
handles.drug_spikes=[];
% initialization for event averages
rate_vect_bis=[];%will contain the rate associated to each event (see method section for the definition of the computation of the rate)
average_length = [];
average_no_spk=[];
required_mean_nb_spikes_per_rendition = 1;
for k=1:length(nodrug_index_events_r_merged)
    t_inf=nodrug_time_events_l_merged(k);
    t_sup=nodrug_time_events_r_merged(k);
    for j=1:length(handles.cell)
        for i=1:handles.noMotifs
            spikeIndex=find(handles.spikes{j,i}>=t_inf & handles.spikes{j,i}<=t_sup);
            test = handles.nodrugIndex == i;
            if sum(test)==0
            handles.nodrug_spikes{j,i}=[-100];
            handles.drug_spikes{j,i}=handles.spikes{j,i}(spikeIndex);
            else
            handles.nodrug_spikes{j,i}=handles.spikes{j,i}(spikeIndex);
            handles.drug_spikes{j,i}=[-100];
            end
        end
    end

[maxi index_max]=max(smooth_psth_nodrug(nodrug_index_events_l_merged(k):nodrug_index_events_r_merged(k)));
time_max=handles.bin(nodrug_index_events_l_merged(k)+index_max-1);
m_left=smooth_psth_nodrug(nodrug_index_events_l_merged(k));
m_right=smooth_psth_nodrug(nodrug_index_events_r_merged(k));
motifTemplate_start=handles.motifTemplate(1);
motifTemplate_end=handles.motifTemplate(end);
t_min_left=nodrug_time_events_l_merged(k);
t_min_right=nodrug_time_events_r_merged(k);
num=k;
num_str=num2str(num);
zero_str=num2str(0);
spk=handles.nodrug_spikes;
spk_bis=handles.drug_spikes;
if num<10
    num_str=strcat(zero_str,num_str);
end
nb_of_max = t_min_left<time_max_nodrug  & time_max_nodrug<t_min_right;
if sum(nb_of_max)==1
    status=1;
elseif sum(nb_of_max)==2 & maxi>20
    status=2; %will identify the double bump events
else
    status=3;%should correspond to low rate events
end
    
seq=num2str(handles.Seq);
nb_nodrug_renditions=num2str(length(handles.nodrugIndex));
file_name_nodrug=strcat('nodrug_',handles.name_cell,'_syll',seq,'_event',num_str);
nodrugIndex=handles.nodrugIndex;
% if num~=(length(time_min_nodrug)-1)
%     plant
% else
save(file_name_nodrug,'spk','spk_bis','nb_nodrug_renditions','t_min_left','t_min_right'...
    ,'time_max','m_left','m_right','maxi','motifTemplate_start','motifTemplate_end'...
    ,'status','nodrugIndex')
%end

%computation of the length of an event and the mean number of spikes per rendition

nb_of_spk = []; %for each event will contain the number of spikes in each rendition
Tend_minus_T1 = [];%for each event will contain the difference between the timing of the last and the first spike in each rendition

    for j=1:length(spk)
        if spk{j}~=[-100] %if for a rendition j, spk{j} = [-100] it means that it corresponds to a rendition during which drug was infused in LMAN
            nb_of_spk = cat(2,nb_of_spk,length(spk{j}));
            if isempty(spk{j})==0
                Tend_minus_T1 = cat(2,Tend_minus_T1,spk{j}(end)-spk{j}(1));
            else
                Tend_minus_T1 = cat(2,Tend_minus_T1,0);
            end
        end
    end
    
    
%-----------
%computation of the rate
    nb_of_rendition=0;
    for j=1:length(spk)
        if spk{j}~=[-100] %loop to deal with the song renditions that are made when drug has been injected
%into LMAN. For those renditions j, the variable spk{j} usually containing
%the timing of the spike is replaced by [-100];
            nb_of_rendition = nb_of_rendition + 1;
        end
    end
    
    mean_duration=[];
    isi_event=[];


    for j=1:length(spk)
        if length(spk{j})>=1 & spk{j}~=[-100]
            if length(spk{j})>=2
               mean_duration = cat(2,mean_duration,spk{j}(end)-spk{j}(1));
            end
        end
    end


    for j=1:length(spk)
        if length(spk{j})>=2
            isi_inside_trial=[];
            for l=2:length(spk{j})
                isi_inside_trial=cat(2,isi_inside_trial,spk{j}(l)-spk{j}(l-1));
            end
                isi_event=cat(2,isi_event,mean(isi_inside_trial));
        elseif length(spk{j})==1 & spk{j}~=[-100]
                isi_event=cat(2,isi_event,mean_duration);
        elseif isempty(spk{j})==1
            isi_event=cat(2,isi_event,inf);
        elseif spk{j}==[-100]
        end
    end
    rate_vect=1./isi_event;
    rate=mean(rate_vect);


% adding all the features of the studied event to the vectors that contain
% statistics for the whole cell. Do not consider events for which the mean
% number of spikes per rendition is below:
% "required_mean_nb_spikes_per_rendition"

    if mean(nb_of_spk) >= required_mean_nb_spikes_per_rendition
        rate_vect_bis = cat(2,rate_vect_bis,rate);
        average_length = cat(2,average_length,mean(Tend_minus_T1));
        average_no_spk=cat(2,average_no_spk,mean(nb_of_spk));
    end
end

handles.event_firing_rate = rate_vect_bis;
handles.event_average_length = average_length;
handles.event_average_spk = average_no_spk;




%{code here}




num_both=0;
for k=1:min(length(nodrug_index_events_r_merged),length(drug_index_events_r_merged))
    vect_diff_r = abs(drug_time_events_r_merged - nodrug_time_events_r_merged(k));
    [min_r index_r] = min(vect_diff_r);
    vect_diff_l=abs(drug_time_events_l_merged - nodrug_time_events_l_merged(k));
    [min_l index_l] = min(vect_diff_l);
    if min_r< 0.015 & min_l<0.015 & index_r==index_l
        num_both=num_both+1;
        num_str=num2str(num_both);
        seq=num2str(handles.Seq);
        file_name_both = strcat('both_',handles.name_cell,'_syll',seq,'_event',num_str);
        t_inf_nodrug=nodrug_time_events_l_merged(k);
    t_sup_nodrug=nodrug_time_events_r_merged(k);
    t_inf_drug=drug_time_events_l_merged(k);
    t_sup_drug=drug_time_events_r_merged(k);
    for j=1:length(handles.cell)
        for i=1:handles.noMotifs
            spikeIndex_drug=find(handles.spikes{j,i}>=t_inf_drug & handles.spikes{j,i}<=t_sup_drug);
            spikeIndex_nodrug=find(handles.spikes{j,i}>=t_inf_nodrug & handles.spikes{j,i}<=t_sup_nodrug);
            test = handles.nodrugIndex == i;
            if sum(test)==0
            handles.nodrug_spikes{j,i}=[-100];
            handles.drug_spikes{j,i}=handles.spikes{j,i}(spikeIndex_drug);
            else
            handles.nodrug_spikes{j,i}=handles.spikes{j,i}(spikeIndex_nodrug);
            handles.drug_spikes{j,i}=[-100];
            end
        end
    end
    

spk_nodrug=handles.nodrug_spikes;
nodrugIndex=handles.nodrugIndex;
[maxi_nodrug index_max_nodrug]=max(smooth_psth_nodrug(nodrug_index_events_l_merged(k):nodrug_index_events_r_merged(k)));
time_max_nodrug=handles.bin(nodrug_index_events_l_merged(k)+index_max_nodrug-1);
m_left_nodrug=smooth_psth_nodrug(nodrug_index_events_l_merged(k));
m_right_nodrug=smooth_psth_nodrug(nodrug_index_events_r_merged(k));
motifTemplate_start=handles.motifTemplate(1);
motifTemplate_end=handles.motifTemplate(end);
t_min_left_nodrug=nodrug_time_events_l_merged(k);
t_min_right_nodrug=nodrug_time_events_r_merged(k);
nb_of_max_nodrug = t_min_left_nodrug<time_max_nodrug  & time_max_nodrug<t_min_right_nodrug;
if sum(nb_of_max_nodrug)==1
    status_nodrug=1;
elseif sum(nb_of_max_nodrug)==2 & maxi_nodrug>20
    status_nodrug=2; %will identify the double bump events
else
    status_nodrug=3;%should correspond to low rate events
end

spk_drug=handles.drug_spikes;
drugIndex=handles.drugIndex;
[maxi_drug index_max_drug]=max(smooth_psth_drug(drug_index_events_l_merged(k):drug_index_events_r_merged(k)));
time_max_drug=handles.bin(drug_index_events_l_merged(k)+index_max_drug-1);
m_left_drug=smooth_psth_drug(drug_index_events_l_merged(k));
m_right_drug=smooth_psth_drug(drug_index_events_r_merged(k));
t_min_left_drug=drug_time_events_l_merged(k);
t_min_right_drug=drug_time_events_r_merged(k);
nb_of_max_drug = t_min_left_drug<time_max_drug  & time_max_drug<t_min_right_drug;
spk=handles.spikes;
if sum(nb_of_max_drug)==1
    status_drug=1;
elseif sum(nb_of_max_drug)==2 & maxi_drug>20
    status_drug=2; %will identify the double bump events
else
    status_drug=3;%should correspond to low rate events
end
        
        save(file_name_both,'spk','spk_drug','spk_nodrug','t_min_left_drug','t_min_right_drug'...
    ,'time_max_drug','m_left_drug','m_right_drug','maxi_drug','t_min_left_nodrug','t_min_right_nodrug'...
    ,'time_max_nodrug','m_left_nodrug','m_right_nodrug','maxi_nodrug','motifTemplate_start','motifTemplate_end'...
    ,'status_drug','status_nodrug','drugIndex','nodrugIndex')


    end
end

%-> Display events
handles.treshold=str2double(get(hObject,'String'));
    subplot(handles.raster);
    for i=1:length(nodrug_time_events_l_merged)
        if isempty(nodrug_boundaries)==0
        for j=1:length(nodrug_boundaries)/2
    fill([nodrug_time_events_l_merged(i) nodrug_time_events_r_merged(i) ...
        nodrug_time_events_r_merged(i) nodrug_time_events_l_merged(i) ]...
        ,[nodrug_boundaries(2*j-1)-1 nodrug_boundaries(2*j-1)-1 nodrug_boundaries(2*j) nodrug_boundaries(2*j)],'g')
    hold on
        end
        end
    end
    for i=1:length(drug_time_events_l_merged)
        if isempty(drug_boundaries)==0
        for j=1:length(drug_boundaries)/2
    fill([drug_time_events_l_merged(i) drug_time_events_r_merged(i) ...
        drug_time_events_r_merged(i) drug_time_events_l_merged(i) ]...
        ,[drug_boundaries(2*j-1)-1 drug_boundaries(2*j-1)-1 drug_boundaries(2*j) drug_boundaries(2*j)],'m')
    hold on
        end
        end
    end
    alpha(0.1);

subplot(handles.spikeStats);
         cla;

                 % plot(handles.bin,handles.psth_nodrug,'r');
         %plot(handles.bin,handles.psth_nodrug,handles.bin,handles.psth_dru
         %g,'r');
         plot(handles.bin,abs(smooth_psth_nodrug),'b');
         hold on
%          for i=1:length(index_max_drug)
%              plot([time_max_drug(i),time_max_drug(i)],[0,400],'g')
%              hold on
%          end
%          for i=1:length(index_min_drug)
%              plot([time_min_drug(i),time_min_drug(i)],[0,400],'y')
%              hold on
%          end


         for i=1:length(nodrug_time_events_l)
             plot([nodrug_time_events_l(i),nodrug_time_events_l(i)],[0,400],'b')
             hold on
         end
         for i=1:length(nodrug_time_events_r)
             ve=nodrug_time_events_l==nodrug_time_events_r(i);
             if ve==0
             plot([nodrug_time_events_r(i),nodrug_time_events_r(i)],[0,400],'r')
             hold on
             else
                 plot([nodrug_time_events_r(i),nodrug_time_events_r(i)],[0,400],'g')
             hold on
             end
                 
         end
       
         
%          plot(handles.bin,handles.treshold*ones(length(handles.bin)))
%         xlim([(handles.motifTemplate(1)-0.1) (handles.motifTemplate(end)+0.1)]);
% set(handles.to_t,'String',num2str(handles.to_time));
end
guidata(hObject, handles);


%--- Executes during object creation, after setting all properties.
function treshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to to (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


