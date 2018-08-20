function varargout = SylClus3(varargin)
% SYLCLUS3 MATLAB code for SylClus3.fig
%      SYLCLUS3, by itself, creates a new SYLCLUS3 or raises the existing
%      singleton*.
%
%      H = SYLCLUS3 returns the handle to a new SYLCLUS3 or the handle to
%      the existing singleton*.
%
%      SYLCLUS3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SYLCLUS3.M with the given input arguments.
%
%      SYLCLUS3('Property','Value',...) creates a new SYLCLUS3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SylClus3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SylClus3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SylClus3

% Last Modified by GUIDE v2.5 18-Nov-2015 11:06:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SylClus3_OpeningFcn, ...
                   'gui_OutputFcn',  @SylClus3_OutputFcn, ...
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

function SylClus3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SylClus3 (see VARARGIN)

% PASTE IN HANDLES
handles.audio = varargin{1}.audio;
handles.sylStart = varargin{1}.sylStart;
handles.sylEnd = varargin{1}.sylEnd;
handles.sylType = varargin{1}.sylType;
handles.IndxKey = varargin{1}.IndxKey;
handles.sylvalue = varargin{1}.sylvalue;
handles.sylPrevalue = varargin{1}.sylPrevalue; %
handles.sylPostvalue = varargin{1}.sylPostvalue; %

handles.syltempvalueNot = varargin{1}.sylvalueNot; %
handles.syltempPreNot = varargin{1}.sylPrevalueNot; %
handles.syltempPostNot = varargin{1}.sylPostvalueNot; %

% PARAMETER FOR THE DISPLAY
%Generate indicesf
[handles.syltempBin, handles.syltempIndx] = createIndices(handles);

handles.syltempStart=handles.sylStart(handles.syltempBin);
handles.syltempEnd=handles.sylEnd(handles.syltempBin);
handles.syltempType=handles.sylType(handles.syltempBin);
handles.syltempIndxKey=handles.IndxKey(handles.syltempBin);

%CLEAR PLOT
axes(handles.axes_syl1); cla
axes(handles.axes_syl2); cla
axes(handles.axes_syl3); cla
axes(handles.axes_syl4); cla

% SPECTRO INTERPOL
SpecSyl=cell(length(handles.syltempStart),1);
FeatSyl=cell(length(handles.syltempStart),1);
entSyl = [];
margin= floor(0.01*44150);
time=zeros(length(SpecSyl),1);

startmargin=floor(handles.syltempStart-margin);
startmargin(startmargin<=0)=1;
endmargin=floor(handles.syltempEnd+margin);

%CALCUL SPECTRO FOR EACH FILE
for i=1:length(handles.syltempStart)
    %Generate snippet of audio corresponding to the syllable
    try
        SpecSyl{i}=handles.audio{handles.syltempIndxKey(i)}(startmargin(i):endmargin(i));
    catch 
        SpecSyl{i}=handles.audio{handles.syltempIndxKey(i)}(startmargin(i):end);
    end
    
    %Create spectrogram for syllable
    [~,~,~,SpecSyl{i}] = spectrogram((SpecSyl{i}/(sqrt(mean(SpecSyl{i}.^2)))),220,220-44,512,44150);
    strt = 4; stp = 94; %These bins correspond to 300-8000Hz
    SpecSyl{i} = abs(SpecSyl{i}(strt:stp,:));
    time(i)=size(SpecSyl{i},2); %find duration of syllable (estimated by the length of the spectrogram)
    
    %Calculate the features of the margin-less snippet
    FeatSyl{i} = koenigSpectral(handles.audio{handles.syltempIndxKey(i)}(floor(handles.syltempStart(i)):ceil(handles.syltempEnd(i))), 44150);
    entSyl(i) = mean(FeatSyl{i}.Entropy);
%     pitchSyl(i) = mean(FeatSyl{i}.Pitch);
    pitchSyl(i) = mean(FeatSyl{i}.Pitch_chose);  
    fmSyl(i) = mean(FeatSyl{i}.FM);
    amSyl(i) = mean(FeatSyl{i}.AM);
    ampSyl(i) = mean(FeatSyl{i}.amplitude);
end
handles.FeatSyl = FeatSyl;
handles.entSyl = entSyl;
handles.pitchSyl = pitchSyl;
handles.fmSyl = fmSyl;
handles.amSyl = amSyl;
handles.ampSyl = ampSyl;
duration=(handles.syltempEnd-handles.syltempStart)/44150;

%PCA decomp
prinMat = [entSyl', pitchSyl', fmSyl', amSyl', ampSyl', duration'];
handles = runPCA(handles, prinMat);

% [Z, handles.mu, handles.sigma] = zscore(prinMat);
% [coeff, score, latent] = princomp(Z);
% handles.pcaCoeff = coeff;
% handles.PC1 = score(:,1)';
% handles.PC2 = score(:,2)';
% handles.PC3 = score(:,3)';

% INTERPOL SPECTRO
freq=size(SpecSyl{1},1);
tim=floor(mean(time));

SpecSylinter=zeros(size(SpecSyl{1},1),tim ,length(SpecSyl));
for i=1:length(SpecSyl)
    for j=1:freq
        SpecSylinter(j,:,i)=interp1(SpecSyl{i}(j,:),linspace(1,size(SpecSyl{i},2),tim));
    end
end
handles.SpecSylinter=SpecSylinter;

% CRITERION TO DISPLAY
% crit_bin=(handles.syltempEnd-handles.syltempStart)>0;
MeanSpec=log10(mean(SpecSylinter,3));

axes(handles.axes_syl1)
imagesc(MeanSpec)
axis xy


% handles.limits=sort([get(handles.slider_bound1,'Value') get(handles.slider_bound2,'Value') get(handles.slider_bound3,'Value')]);
% duration=(handles.syltempEnd-handles.syltempStart)/44150;

% SET SLIDER AND EDIT BOX
%%%%%%%%%%%%%%%%%%%%
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
set(handles.slider_split1, 'SliderStep', [0.05 1]) 

set(handles.slider_split2,'Value',0)
set(handles.slider_split2,'Min',0)
set(handles.slider_split2,'Max',tim)
set(handles.slider_split2, 'SliderStep', [0.05 1])

set(handles.slider_split3,'Value',0)
set(handles.slider_split3,'Min',0)
set(handles.slider_split3,'Max',tim)
set(handles.slider_split3, 'SliderStep', [0.05 1]) 

set(handles.slider_split4,'Value',0)
set(handles.slider_split4,'Min',0)
set(handles.slider_split4,'Max',tim)
set(handles.slider_split4, 'SliderStep', [0.05 1]) 

% DISPLAY HISTO
handles=refresh_hist(handles);

% Choose default command line output for SylClus3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SylClus3 wait for user response (see UIRESUME)
uiwait(handles.figure1);

function varargout = SylClus3_OutputFcn(hObject, eventdata, handles)

%Format data to be sent to the calling GUI
output_sylclus.sylStart=handles.sylStart ;
output_sylclus.sylEnd=handles.sylEnd ;
output_sylclus.sylType=handles.sylType ;
output_sylclus.IndxKey=handles.IndxKey ;

varargout{1}=output_sylclus;
delete(handles.figure1);

function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Independent functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles=refresh_hist(handles)    
%Clear the current figure
axes(handles.axes_hist); cla

%Retrieve User selection for X axis
selX = get(handles.popup_featureX,'Value');
if selX == 1
    featX = handles.entSyl;
elseif selX == 2
    featX = handles.pitchSyl;
elseif selX == 3
    featX = handles.fmSyl;
elseif selX == 4
    featX = handles.amSyl;
elseif selX == 5
    featX = handles.ampSyl;
elseif selX == 6
    featX = handles.PC1;
elseif selX == 7
    featX = handles.PC2;
elseif selX == 8
    featX = handles.PC3;
elseif selX == 9
    featX = (handles.syltempEnd-handles.syltempStart)/44150;
end

%Retrieve User selection for Y axis
selY = get(handles.popup_featureY,'Value');
if selY == 1
    featY = handles.entSyl;
elseif selY == 2
    featY = handles.pitchSyl;
elseif selY == 3
    featY = handles.fmSyl;
elseif selY == 4
    featY = handles.amSyl;
elseif selY == 5
    featY = handles.ampSyl;
elseif selY == 6
    featY = handles.PC1;
elseif selY == 7
    featY = handles.PC2;
elseif selY == 8
    featY = handles.PC3;
elseif selY == 9
    featY = (handles.syltempEnd-handles.syltempStart)/44150;
end

%Do the histogram
if length(featX) ~= 1
    [handles.edgesX2,handles.edgesY2,handles.N,handles.histImgHandle] = ndhist(featX, featY, 'bins', 10, 'filter', 'max');
else
    scatter(featX, featY, 'xk');
end

box off
set(gca, 'TickDir', 'out')
xlabel('Feature X', 'FontSize', 20)
ylabel('Feature Y', 'FontSize', 20)
set(gca, 'LineWidth', 3, 'FontSize', 14)
legend('boxoff')

function [syltempBin, syltempIndx] = createIndices(handles)
%Generate a binary index and temporary sylType array for the selected syllables.
 if ~isnan(handles.sylvalue)
    syltarBin = handles.sylType==handles.sylvalue;
else
    syltarBin = true(size(handles.sylType)); % wildcard
 end

if handles.syltempvalueNot
    syltarBin = ~syltarBin;
end
 

if ~isnan(handles.sylPrevalue)
    sylpreBin = handles.sylType==handles.sylPrevalue;
else
    sylpreBin = true(size(syltarBin)); % wildcard
end

if handles.syltempPreNot
    sylpreBin = ~sylpreBin;
end

if ~isnan(handles.sylPostvalue)
    sylpostBin = handles.sylType==handles.sylPostvalue;
else
    sylpostBin = true(size(syltarBin)); % wildcard
end

if handles.syltempPostNot
    sylpostBin = ~sylpostBin;
end

ind = 2:(length(syltarBin)-1);
syltempBin = syltarBin(ind) & sylpreBin(ind-1) & sylpostBin(ind+1);
syltempBin = [false; syltempBin; false];
if isnan(handles.sylPrevalue) && syltarBin(1) == true
    syltempBin(1) = true;
end
if isnan(handles.sylPostvalue) && syltarBin(end) == true
    syltempBin(end) = true;
end

%Generate a pointer for each syllable back to the main sylType array
syltempIndx = find(syltempBin == true);

function ROIposition = drawROI(handles)
%Draw ROI using mouse and collect outline points. Must be on separate figure because of some weird MATLAB qwerk
f1 = figure; 
% sel = get(handles.popup_featureY,'Value');
% if sel == 1
%     feat = handles.entSyl;
% elseif sel == 2
%     feat = handles.pitchSyl;
% elseif sel == 3
%     feat = handles.fmSyl;
% elseif sel == 4
%     feat = handles.amSyl;
% elseif sel == 5
%     feat = handles.ampSyl;
% end
% 
% if length(feat) ~= 1
%     ndhist((handles.syltempEnd-handles.syltempStart)/44150, feat, 'bins', 10, 'filter', 'max');
% else
%     scatter((handles.syltempEnd-handles.syltempStart)/44150, feat, 'xk');
% end

%Retrieve User selection for X axis
selX = get(handles.popup_featureX,'Value');
if selX == 1
    featX = handles.entSyl;
elseif selX == 2
    featX = handles.pitchSyl;
elseif selX == 3
    featX = handles.fmSyl;
elseif selX == 4
    featX = handles.amSyl;
elseif selX == 5
    featX = handles.ampSyl;
elseif selX == 6
    featX = handles.PC1;
elseif selX == 7
    featX = handles.PC2;
elseif selX == 8
    featX = handles.PC3;
elseif selX == 9
    featX = (handles.syltempEnd-handles.syltempStart)/44150;
end

%Retrieve User selection for Y axis
selY = get(handles.popup_featureY,'Value');
if selY == 1
    featY = handles.entSyl;
elseif selY == 2
    featY = handles.pitchSyl;
elseif selY == 3
    featY = handles.fmSyl;
elseif selY == 4
    featY = handles.amSyl;
elseif selY == 5
    featY = handles.ampSyl;
elseif selY == 6
    featY = handles.PC1;
elseif selY == 7
    featY = handles.PC2;
elseif selY == 8
    featY = handles.PC3;
elseif selY == 9
    featY = (handles.syltempEnd-handles.syltempStart)/44150;
end

%Do the histogram
if length(featX) ~= 1
    ndhist(featX, featY, 'bins', 10, 'filter', 'max');
else
    scatter(featX, featY, 'xk');
end

ROI = imfreehand(gca);
ROIposition = getPosition(ROI);
close(f1)

function selInd = getROIsubset(handles, mask)
% Use the mask to generate an index of syllable that fall within the ROI

%Retrieve User selection for X axis
selX = get(handles.popup_featureX,'Value');
if selX == 1
    featX = handles.entSyl;
elseif selX == 2
    featX = handles.pitchSyl;
elseif selX == 3
    featX = handles.fmSyl;
elseif selX == 4
    featX = handles.amSyl;
elseif selX == 5
    featX = handles.ampSyl;
elseif selX == 6
    featX = handles.PC1;
elseif selX == 7
    featX = handles.PC2;
elseif selX == 8
    featX = handles.PC3;
elseif selX == 9
    featX = (handles.syltempEnd-handles.syltempStart)/44150;
end

%Retrieve User selection for Y axis
selY = get(handles.popup_featureY,'Value');
if selY == 1
    featY = handles.entSyl;
elseif selY == 2
    featY = handles.pitchSyl;
elseif selY == 3
    featY = handles.fmSyl;
elseif selY == 4
    featY = handles.amSyl;
elseif selY == 5
    featY = handles.ampSyl;
elseif selY == 6
    featY = handles.PC1;
elseif selY == 7
    featY = handles.PC2;
elseif selY == 8
    featY = handles.PC3;
elseif selY == 9
    featY = (handles.syltempEnd-handles.syltempStart)/44150;
end

%Generate index of points that fall within the mask
for i =1:length(featX)
    [~, pntsx(i)] = min(abs((handles.edgesX2-featX(i))));
    [~, pntsy(i)] = min(abs((handles.edgesY2-featY(i))));
    selInd(i) = mask(pntsy(i),pntsx(i));
end

function handles = resetROIs(handles)
%Reload the orginal value into the edit boxes
set(handles.edit_rename1,'String',handles.sylvalue)
set(handles.edit_rename2,'String',handles.sylvalue)
set(handles.edit_rename3,'String',handles.sylvalue)
set(handles.edit_rename4,'String',handles.sylvalue)

%Clear off the meanSyl axes
axes(handles.axes_syl2); cla
axes(handles.axes_syl3); cla
axes(handles.axes_syl4); cla

%Delete the masks ROI2position
if isfield(handles, 'mask1')
    handles = rmfield(handles, {'mask1', 'ROI1position', 'ROI1Ind'});
end
if isfield(handles, 'mask2')
    handles = rmfield(handles, {'mask2', 'ROI2position', 'ROI2Ind'});
end
if isfield(handles, 'mask3')
    handles = rmfield(handles, {'mask3', 'ROI3position', 'ROI3Ind'});
end

function handles = runPCA(handles, prinMat)
%Run the PCA decomposition on the current feature set
%prinMat = [entSyl', pitchSyl', fmSyl', amSyl', ampSyl', duration'];

%Zscore the feature set
[Z, handles.mu, handles.sigma] = zscore(prinMat);

%Do PCA on the z-scored data
[coeff, score, latent] = princomp(Z);

%Copy out the useful info
handles.pcaCoeff = coeff;
handles.latent = latent;
handles.PC1 = score(:,1)';
handles.PC2 = score(:,2)';
handles.PC3 = score(:,3)';

%Display variance explanation plot
axes(handles.axes_pcaExp); cla
varExp = cumsum(latent)./sum(latent);
plot(1:numel(varExp), varExp, '--sr', 'LineWidth', 1);
xlabel('PC Number', 'FontSize', 10); xlim([0.5, numel(varExp)+0.5])
ylabel('%', 'FontSize', 10); ylim([0,1]);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 1:numel(varExp), 'LineWidth', 2, 'FontSize', 10)

%Display Factor loadings
axes(handles.axes_pcaLoad); cla
bar(handles.pcaCoeff')
xlabel('PC Number', 'FontSize', 10); xlim([0.5, numel(varExp)+0.5])
ylabel('Loading', 'FontSize', 10);
set(gca, 'Box', 'off', 'TickDir', 'out','XTick', 1:numel(varExp), 'LineWidth', 2, 'FontSize', 10)
legend('Ent', 'Pitch','FM', 'AM', 'Amp', 'Dur')
legend([650,5,10,5], 'Orientation', 'horizontal'); legend('boxoff')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histogram display tools
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_pcaAgain_Callback(hObject, eventdata, handles)
%Run the PCA decomposition on the current feature set

%Generate indices
[handles.syltempBin, handles.syltempIndx] = createIndices(handles);

handles.syltempStart=handles.sylStart(handles.syltempBin);
handles.syltempEnd=handles.sylEnd(handles.syltempBin);
handles.syltempType=handles.sylType(handles.syltempBin);
handles.syltempIndxKey=handles.IndxKey(handles.syltempBin);

%Reset ROI stuff
handles = resetROIs(handles);

%Run PCA
%Make the features matrix
duration=(handles.syltempEnd-handles.syltempStart)/44150;
prinMat = [handles.entSyl', handles.pitchSyl', handles.fmSyl', handles.amSyl', handles.ampSyl', duration'];

handles = runPCA(handles, prinMat);

%Redisplay the sorting image
handles=refresh_hist(handles);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROI Selection & Reassignment Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function push_refresh_Callback(hObject, eventdata, handles)
%Generate indices
[handles.syltempBin, handles.syltempIndx] = createIndices(handles);

handles.syltempStart=handles.sylStart(handles.syltempBin);
handles.syltempEnd=handles.sylEnd(handles.syltempBin);
handles.syltempType=handles.sylType(handles.syltempBin);
handles.syltempIndxKey=handles.IndxKey(handles.syltempBin);

%Update those sylTypes that have been defined by the ROIs
%If there is a mask for ROI1, remap sylTypes as specified
if isfield(handles, 'mask1') && (str2double(get(handles.edit_rename2,'String')) ~= handles.sylvalue)
    % Use the mask to generate an index of syllable that fall within the ROI
    ROI1Ind = getROIsubset(handles, handles.mask1);
    
    %Relabel the types
    handles.sylType(handles.syltempIndx(ROI1Ind)) = str2double(get(handles.edit_rename2,'String'));
else
    ROI1Ind = false(size(handles.syltempIndx))';
end

%If there is a mask for ROI2, remap sylTypes as specified
if isfield(handles, 'mask2') && (str2double(get(handles.edit_rename3,'String')) ~= handles.sylvalue)
    % Use the mask to generate an index of syllable that fall within the ROI
    ROI2Ind = getROIsubset(handles, handles.mask2);
    
    %Relabel the types
    handles.sylType(handles.syltempIndx(ROI2Ind)) = str2double(get(handles.edit_rename3,'String'));
else
    ROI2Ind = false(size(handles.syltempIndx))';
end

%If there is a mask for ROI3, remap sylTypes as specified
if isfield(handles, 'mask3') && (str2double(get(handles.edit_rename4,'String')) ~= handles.sylvalue)
    % Use the mask to generate an index of syllable that fall within the ROI
    ROI3Ind = getROIsubset(handles, handles.mask3);
    
    %Relabel the types
    handles.sylType(handles.syltempIndx(ROI3Ind)) = str2double(get(handles.edit_rename4,'String'));
else
    ROI3Ind = false(size(handles.syltempIndx))';
end

%Remove the corresponding data from the feature arrays
sumROI = (ROI1Ind | ROI2Ind) | ROI3Ind;
handles.FeatSyl(sumROI) = [];
handles.entSyl(sumROI) = [];
handles.pitchSyl(sumROI) = [];
handles.fmSyl(sumROI) = [];
handles.amSyl(sumROI) = [];
handles.ampSyl(sumROI) = [];
handles.PC1(sumROI) = [];
handles.PC2(sumROI) = [];
handles.PC3(sumROI) = [];

%And from the spectrogram file
handles.SpecSylinter(:,:,sumROI) = [];

%Reset ROI stuff
handles = resetROIs(handles);

%Generate indices
[handles.syltempBin, handles.syltempIndx] = createIndices(handles);

handles.syltempStart=handles.sylStart(handles.syltempBin);
handles.syltempEnd=handles.sylEnd(handles.syltempBin);
handles.syltempType=handles.sylType(handles.syltempBin);
handles.syltempIndxKey=handles.IndxKey(handles.syltempBin);

%Redisplay the sorting image
handles=refresh_hist(handles);

%Update the Ungrouped mean syl plot
MeanSpec=log10(mean(handles.SpecSylinter(:,:,:),3));
axes(handles.axes_syl1)
imagesc(MeanSpec)
axis xy

guidata(hObject, handles);

function push_ROI1_Callback(hObject, eventdata, handles)
%If present, delete the previous ROI drawing on the axis
if isfield(handles, 'ROI1plot')
    try
        delete(handles.ROI1plot)
    catch
        handles = rmfield(handles, 'ROI1plot');
    end
end

%Draw ROI and gather points defining the outline
[handles.ROI1position] = drawROI(handles);

%Translate the position of the ROI into image coordinates
for i = 1:size(handles.ROI1position,1)
    [~, indx(i)] = min(abs((handles.edgesX2-handles.ROI1position(i,1))));
    [~, indy(i)] = min(abs((handles.edgesY2-handles.ROI1position(i,2))));
end

%Generate image mask from transformed points 
mask = poly2mask(indx, indy, length(handles.edgesY2), length(handles.edgesX2));

%Check that this mask doesn't overlap with any other mask; if it does, bounce it
m2OK = false;
if isfield(handles, 'mask2') %Check mask 1 for overlap
    if ~any(handles.mask2(:) & mask(:))
        m2OK = true;
    end
else
    m2OK = true;
end

m3OK = false;
if isfield(handles, 'mask3') %Check mask 2 overlap
    if ~any(handles.mask3(:) & mask(:))
        m3OK = true;
    end
else
    m3OK = true;
end

if m2OK && m3OK
    %No overlap found, so generate mask and continue
    handles.mask1 = mask;
else
    %Overlaps found, so trash mask and clear out old data
    if isfield(handles, 'mask1')
        handles = rmfield(handles,'mask1');
    end
    axes(handles.axes_syl2); cla;
    warndlg('ROI overlaps with previously defined ROI. Try again.')
    return
end

%Plot the new ROI on the main axis
axes(handles.axes_hist)
hold on
handles.ROI1plot = plot(handles.ROI1position(:,1),handles.ROI1position(:,2), 'r');
hold off

% Use the mask to generate an index of syllable that fall within the ROI
handles.ROI1Ind = getROIsubset(handles, handles.mask1);

% Update the composite syllable spectrogram
handles.MeanSpec(:,:,2)=log10(mean(handles.SpecSylinter(:, :, handles.ROI1Ind), 3));
axes(handles.axes_syl2)
imagesc(handles.MeanSpec(:,:,2))
axis xy

guidata(hObject, handles);

function push_ROI2_Callback(hObject, eventdata, handles)
%If present, delete the previous ROI drawing on the axis
if isfield(handles, 'ROI2plot')
    try
        delete(handles.ROI2plot)
    catch
        handles = rmfield(handles, 'ROI2plot');
    end
end

%Draw ROI and gather points defining the outline
[handles.ROI2position] = drawROI(handles);

%Translate the position of the ROI into image coordinates
for i = 1:size(handles.ROI2position,1)
    [~, indx(i)] = min(abs((handles.edgesX2-handles.ROI2position(i,1))));
    [~, indy(i)] = min(abs((handles.edgesY2-handles.ROI2position(i,2))));
end

%Generate image mask from transformed points 
mask = poly2mask(indx, indy, length(handles.edgesY2), length(handles.edgesX2));

%Check that this mask doesn't overlap with any other mask; if it does, bounce it
m1OK = false;
if isfield(handles, 'mask1') %Check mask 1 for overlap
    if ~any(handles.mask1(:) & mask(:))
        m1OK = true;
    end
else
    m1OK = true;
end

m3OK = false;
if isfield(handles, 'mask3') %Check mask 2 overlap
    if ~any(handles.mask3(:) & mask(:))
        m3OK = true;
    end
else
    m3OK = true;
end

if m1OK && m3OK
    %No overlap found, so generate mask and continue
    handles.mask2 = mask;
else
    %Overlaps found, so trash mask and clear out old data
    if isfield(handles, 'mask2')
        handles = rmfield(handles,'mask2');
    end
        axes(handles.axes_syl3); cla;
    warndlg('ROI overlaps with previously defined ROI. Try again.')
    return
end

%Plot the new ROI on the main axis
axes(handles.axes_hist)
hold on
handles.ROI2plot = plot(handles.ROI2position(:,1),handles.ROI2position(:,2), 'g');
hold off

% Use the mask to generate an index of syllable that fall within the ROI
handles.ROI2Ind = getROIsubset(handles, handles.mask2);

% Update the composite syllable spectrogram
handles.MeanSpec(:,:,3)=log10(mean(handles.SpecSylinter(:, :, handles.ROI2Ind), 3));
axes(handles.axes_syl3)
imagesc(handles.MeanSpec(:,:,3))
axis xy

guidata(hObject, handles);

function push_ROI3_Callback(hObject, eventdata, handles)
%If present, delete the previous ROI drawing on the axis
if isfield(handles, 'ROI3plot')
    try
        delete(handles.ROI3plot)
    catch
        handles = rmfield(handles, 'ROI3plot');
    end
end

%Draw ROI and gather points defining the outline
[handles.ROI3position] = drawROI(handles);

%Translate the position of the ROI into image coordinates
for i = 1:size(handles.ROI3position,1)
    [~, indx(i)] = min(abs((handles.edgesX2-handles.ROI3position(i,1))));
    [~, indy(i)] = min(abs((handles.edgesY2-handles.ROI3position(i,2))));
end

%Generate image mask from transformed points 
mask = poly2mask(indx, indy, length(handles.edgesY2), length(handles.edgesX2));

%Check that this mask doesn't overlap with any other mask; if it does, bounce it
m1OK = false;
if isfield(handles, 'mask1') %Check mask 1 for overlap
    if ~any(handles.mask1(:) & mask(:))
        m1OK = true;
    end
else
    m1OK = true;
end

m2OK = false;
if isfield(handles, 'mask2') %Check mask 2 overlap
    if ~any(handles.mask2(:) & mask(:))
        m2OK = true;
    end
else
    m2OK = true;
end

if m1OK && m2OK
    %No overlap found, so generate mask and continue
    handles.mask3 = mask;
else
    %Overlaps found, so trash mask and clear out old data
    if isfield(handles, 'mask3')
        handles = rmfield(handles,'mask3');
    end
    axes(handles.axes_syl4); cla;
    warndlg('ROI overlaps with previously defined ROI. Try again.')
    return
end

%Plot the new ROI on the main axis
axes(handles.axes_hist)
hold on
handles.ROI3plot = plot(handles.ROI3position(:,1),handles.ROI3position(:,2), 'k');
hold off

% Use the mask to generate an index of syllable that fall within the ROI
handles.ROI3Ind = getROIsubset(handles, handles.mask3);

% Update the composite syllable spectrogram
handles.MeanSpec(:,:,4)=log10(mean(handles.SpecSylinter(:, :, handles.ROI3Ind), 3));
axes(handles.axes_syl4)
imagesc(handles.MeanSpec(:,:,4))
axis xy

guidata(hObject, handles);

function popup_featureY_Callback(hObject, eventdata, handles)
%Reset ROI stuff
handles = resetROIs(handles);

%Refresh the histogram
handles=refresh_hist(handles);

guidata(hObject, handles);

function popup_featureY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_featureX_Callback(hObject, eventdata, handles)
%Reset ROI stuff
handles = resetROIs(handles);

%Refresh the histogram
handles=refresh_hist(handles);

guidata(hObject, handles);

function popup_featureX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function push_split1_Callback(hObject, eventdata, handles)

function push_split2_Callback(hObject, eventdata, handles)
%Generate indices
[handles.syltempBin, handles.syltempIndx] = createIndices(handles);

handles.syltempStart=handles.sylStart(handles.syltempBin);
handles.syltempEnd=handles.sylEnd(handles.syltempBin);
handles.syltempType=handles.sylType(handles.syltempBin);
handles.syltempIndxKey=handles.IndxKey(handles.syltempBin);

%Update only those sylTypes that have been defined by the ROIs
%If there is a mask for ROI1, remap sylTypes as specified
if isfield(handles, 'mask1')
    % Use the mask to generate an index of syllable that fall within the ROI
    ROI1Ind = getROIsubset(handles, handles.mask1);
    localPntr = find(ROI1Ind);
    
    % Calculate user-specified splitting threshold in %
    value=get(handles.slider_split2,'Value')/size(handles.MeanSpec(:,:,2),2);

    % Calculate the duration of ROI syllables in seconds based on annotion times
    duration=(handles.syltempEnd(ROI1Ind) - handles.syltempStart(ROI1Ind))/44150;
    margin=100;
    
    % Do the splitting
    %Loop through each syllable in the ROI
    ind = 0;
    a = handles.syltempIndx(ROI1Ind); %creates global pointer array
    for i = 1:length(a)
        j=a(i)+i-1; %Requires sequencing... do not parallelize
        
        %Calculate the time (in samples) from syllable onset to end of left syllable
        valueaddstart=floor(duration(i)*44150*value+margin);
        
        %Insert new points into the global arrays
        handles.sylStart=horzcat(handles.sylStart(1:j), handles.sylStart(j)+valueaddstart, handles.sylStart(j+1:end));
        handles.sylEnd=horzcat(handles.sylEnd(1:j-1), handles.sylStart(j+1)-floor(2*margin), handles.sylEnd(j:end));
        handles.sylType=vertcat(handles.sylType(1:j-1), str2double(get(handles.edit_split_left2,'String')), str2double(get(handles.edit_split_right2,'String')), handles.sylType(j+1:end));
        handles.IndxKey=vertcat(handles.IndxKey(1:j), handles.IndxKey(j), handles.IndxKey(j+1:end));

        %Update the feature arrays
%         expLocal = localPntr(i) + i - 1;
        expLocal = localPntr(i) + ind;
        ind = ind + 1;
        tempTot = handles.FeatSyl{expLocal};
        snipLength = length(tempTot.Pitch);
        tempL.AM = tempTot.AM(1:floor(snipLength*value)); tempR.AM =  tempTot.AM((floor(snipLength*value)+1):end);
        tempL.FM = tempTot.FM(1:floor(snipLength*value)); tempR.FM =  tempTot.FM((floor(snipLength*value)+1):end);
        tempL.Entropy = tempTot.Entropy(1:floor(snipLength*value)); tempR.Entropy =  tempTot.Entropy((floor(snipLength*value)+1):end);
        tempL.amplitude = tempTot.amplitude(1:floor(snipLength*value)); tempR.amplitude =  tempTot.amplitude((floor(snipLength*value)+1):end);
        tempL.gravity_center = tempTot.gravity_center(1:floor(snipLength*value)); tempR.gravity_center =  tempTot.gravity_center((floor(snipLength*value)+1):end);
        tempL.PitchGoodness = tempTot.PitchGoodness(1:floor(snipLength*value)); tempR.PitchGoodness =  tempTot.PitchGoodness((floor(snipLength*value)+1):end);
        tempL.Pitch = tempTot.Pitch(1:floor(snipLength*value)); tempR.Pitch =  tempTot.Pitch((floor(snipLength*value)+1):end);

        entL = mean(tempL.Entropy); entR = mean(tempR.Entropy);
        pitchL = mean(tempL.Pitch); pitchR = mean(tempR.Pitch);
        fmL = mean(tempL.FM); fmR = mean(tempR.FM);
        amL = mean(tempL.AM); amR = mean(tempR.AM);
        ampL = mean(tempL.amplitude); ampR = mean(tempR.amplitude);
        durL = length(tempL.Entropy)/1000; durR = length(tempR.Entropy)/1000;
        
        %Create PCA projections for the new elements
        fL = [entL, pitchL, fmL, amL, ampL, durL]; fR = [entR, pitchR, fmR, amR, ampR, durR];
        fL_z = (fL-handles.mu)./handles.sigma; fR_z = (fR-handles.mu)./handles.sigma;
        PC1L = fL_z/handles.pcaCoeff(:,1)'; PC1R = fR_z/handles.pcaCoeff(:,1)';
        PC2L = fL_z/handles.pcaCoeff(:,2)'; PC2R = fR_z/handles.pcaCoeff(:,2)';
        PC3L = fL_z/handles.pcaCoeff(:,3)'; PC3R = fR_z/handles.pcaCoeff(:,3)';
        
        %Split the corresponding Spectrogram
        tempSpec = handles.SpecSylinter(:,:,expLocal);
        specDur = size(tempSpec,2);
        specLorig = tempSpec(:,1:(floor(specDur*value)+3)); specRorig = tempSpec(:,(floor(specDur*value)-3):end);
        for k = 1:size(tempSpec,1)
            specL(k,:) = interp1(1:size(specLorig,2), specLorig(k,:), linspace(1, size(specLorig,2), specDur));
            specR(k,:) = interp1(1:size(specRorig,2), specRorig(k,:), linspace(1, size(specRorig,2), specDur));
        end

        if str2double(get(handles.edit_split_left2,'String')) ~= handles.sylvalue
            tempL = []; entL = []; pitchL = []; fmL = []; amL = []; ampL = []; specL = []; PC1L = []; PC2L = []; PC3L = []; 
            ind = ind - 1;
        end

        if str2double(get(handles.edit_split_right2,'String')) ~= handles.sylvalue
            tempR = []; entR = []; pitchR = []; fmR = []; amR = []; ampR = []; specR = []; PC1R = []; PC2R = []; PC3R = []; 
            ind = ind - 1;
        end

        %Do insertions
        handles.FeatSyl = vertcat(handles.FeatSyl(1:(expLocal-1)), tempL, tempR, handles.FeatSyl((expLocal+1):end));
        handles.entSyl = horzcat(handles.entSyl(1:(expLocal-1)), entL, entR, handles.entSyl((expLocal+1):end));
        handles.pitchSyl = horzcat(handles.pitchSyl(1:(expLocal-1)), pitchL, pitchR, handles.pitchSyl((expLocal+1):end));
        handles.fmSyl = horzcat(handles.fmSyl(1:(expLocal-1)), fmL, fmR, handles.fmSyl((expLocal+1):end));
        handles.amSyl = horzcat(handles.amSyl(1:(expLocal-1)), amL, amR, handles.amSyl((expLocal+1):end));
        handles.ampSyl = horzcat(handles.ampSyl(1:(expLocal-1)), ampL, ampR, handles.ampSyl((expLocal+1):end));
        handles.SpecSylinter = cat(3,handles.SpecSylinter(:,:,1:(expLocal-1)), specL, specR, handles.SpecSylinter(:,:,(expLocal+1):end));
        handles.PC1 = horzcat(handles.PC1(1:(expLocal-1)), PC1L, PC1R, handles.PC1((expLocal+1):end));
        handles.PC2 = horzcat(handles.PC2(1:(expLocal-1)), PC2L, PC2R, handles.PC2((expLocal+1):end));
        handles.PC3 = horzcat(handles.PC3(1:(expLocal-1)), PC3L, PC3R, handles.PC3((expLocal+1):end));
    end
end

%Reset ROI stuff
handles = resetROIs(handles);

%Generate indices
[handles.syltempBin, handles.syltempIndx] = createIndices(handles);

handles.syltempStart=handles.sylStart(handles.syltempBin);
handles.syltempEnd=handles.sylEnd(handles.syltempBin);
handles.syltempType=handles.sylType(handles.syltempBin);
handles.syltempIndxKey=handles.IndxKey(handles.syltempBin);

%Redisplay the sorting image
handles=refresh_hist(handles);

%Update the Ungrouped mean syl plot
MeanSpec=log10(mean(handles.SpecSylinter(:,:,:),3));
axes(handles.axes_syl1)
imagesc(MeanSpec)
axis xy

guidata(hObject, handles);

function push_split3_Callback(hObject, eventdata, handles)
%Generate indices
[handles.syltempBin, handles.syltempIndx] = createIndices(handles);

handles.syltempStart=handles.sylStart(handles.syltempBin);
handles.syltempEnd=handles.sylEnd(handles.syltempBin);
handles.syltempType=handles.sylType(handles.syltempBin);
handles.syltempIndxKey=handles.IndxKey(handles.syltempBin);

%Update only those sylTypes that have been defined by the ROIs
%If there is a mask for ROI1, remap sylTypes as specified
if isfield(handles, 'mask2')
    % Use the mask to generate an index of syllable that fall within the ROI
    ROI2Ind = getROIsubset(handles, handles.mask2);
    localPntr = find(ROI2Ind);
    
    % Calculate user-specified splitting threshold in %
    value=get(handles.slider_split3,'Value')/size(handles.MeanSpec(:,:,3),2);

    % Calculate the duration of ROI syllables in seconds based on annotion times
    duration=(handles.syltempEnd(ROI2Ind) - handles.syltempStart(ROI2Ind))/44150;
    margin=100;
    
    % Do the splitting
    %Loop through each syllable in the ROI
    ind = 0;
    a = handles.syltempIndx(ROI2Ind); %creates global pointer array
    for i = 1:length(a)
        j=a(i)+i-1; %Requires sequencing... do not parallelize
        
        %Calculate the time (in samples) from syllable onset to end of left syllable
        valueaddstart=floor(duration(i)*44150*value+margin);
        
        %Insert new points into the global arrays
        handles.sylStart=horzcat(handles.sylStart(1:j), handles.sylStart(j)+valueaddstart, handles.sylStart(j+1:end));
        handles.sylEnd=horzcat(handles.sylEnd(1:j-1), handles.sylStart(j+1)-floor(2*margin), handles.sylEnd(j:end));
        handles.sylType=vertcat(handles.sylType(1:j-1), str2double(get(handles.edit_split_left3,'String')), str2double(get(handles.edit_split_right3,'String')), handles.sylType(j+1:end));
        handles.IndxKey=vertcat(handles.IndxKey(1:j), handles.IndxKey(j), handles.IndxKey(j+1:end));

        %Update the feature arrays
%         expLocal = localPntr(i) + i - 1;
        expLocal = localPntr(i) + ind;
        ind = ind + 1;
        tempTot = handles.FeatSyl{expLocal};
        snipLength = length(tempTot.Pitch);
        tempL.AM = tempTot.AM(1:floor(snipLength*value)); tempR.AM =  tempTot.AM((floor(snipLength*value)+1):end);
        tempL.FM = tempTot.FM(1:floor(snipLength*value)); tempR.FM =  tempTot.FM((floor(snipLength*value)+1):end);
        tempL.Entropy = tempTot.Entropy(1:floor(snipLength*value)); tempR.Entropy =  tempTot.Entropy((floor(snipLength*value)+1):end);
        tempL.amplitude = tempTot.amplitude(1:floor(snipLength*value)); tempR.amplitude =  tempTot.amplitude((floor(snipLength*value)+1):end);
        tempL.gravity_center = tempTot.gravity_center(1:floor(snipLength*value)); tempR.gravity_center =  tempTot.gravity_center((floor(snipLength*value)+1):end);
        tempL.PitchGoodness = tempTot.PitchGoodness(1:floor(snipLength*value)); tempR.PitchGoodness =  tempTot.PitchGoodness((floor(snipLength*value)+1):end);
        tempL.Pitch = tempTot.Pitch(1:floor(snipLength*value)); tempR.Pitch =  tempTot.Pitch((floor(snipLength*value)+1):end);

        entL = mean(tempL.Entropy); entR = mean(tempR.Entropy);
        pitchL = mean(tempL.Pitch); pitchR = mean(tempR.Pitch);
        fmL = mean(tempL.FM); fmR = mean(tempR.FM);
        amL = mean(tempL.AM); amR = mean(tempR.AM);
        ampL = mean(tempL.amplitude); ampR = mean(tempR.amplitude);
        durL = length(tempL.Entropy)/1000; durR = length(tempR.Entropy)/1000;
        
        %Create PCA projections for the new elements
        fL = [entL, pitchL, fmL, amL, ampL, durL]; fR = [entR, pitchR, fmR, amR, ampR, durR];
        fL_z = (fL-handles.mu)./handles.sigma; fR_z = (fR-handles.mu)./handles.sigma;
        PC1L = fL_z/handles.pcaCoeff(:,1)'; PC1R = fR_z/handles.pcaCoeff(:,1)';
        PC2L = fL_z/handles.pcaCoeff(:,2)'; PC2R = fR_z/handles.pcaCoeff(:,2)';
        PC3L = fL_z/handles.pcaCoeff(:,3)'; PC3R = fR_z/handles.pcaCoeff(:,3)';
        
        %Split the corresponding Spectrogram
        tempSpec = handles.SpecSylinter(:,:,expLocal);
        specDur = size(tempSpec,2);
        specLorig = tempSpec(:,1:(floor(specDur*value)+3)); specRorig = tempSpec(:,(floor(specDur*value)-3):end);
        for k = 1:size(tempSpec,1)
            specL(k,:) = interp1(1:size(specLorig,2), specLorig(k,:), linspace(1, size(specLorig,2), specDur));
            specR(k,:) = interp1(1:size(specRorig,2), specRorig(k,:), linspace(1, size(specRorig,2), specDur));
        end

        if str2double(get(handles.edit_split_left3,'String')) ~= handles.sylvalue
            tempL = []; entL = []; pitchL = []; fmL = []; amL = []; ampL = []; specL = []; PC1L = []; PC2L = []; PC3L = []; 
            ind = ind - 1;
        end

        if str2double(get(handles.edit_split_right3,'String')) ~= handles.sylvalue
            tempR = []; entR = []; pitchR = []; fmR = []; amR = []; ampR = []; specR = []; PC1R = []; PC2R = []; PC3R = []; 
            ind = ind - 1;
        end

        %Do insertions
        handles.FeatSyl = vertcat(handles.FeatSyl(1:(expLocal-1)), tempL, tempR, handles.FeatSyl((expLocal+1):end));
        handles.entSyl = horzcat(handles.entSyl(1:(expLocal-1)), entL, entR, handles.entSyl((expLocal+1):end));
        handles.pitchSyl = horzcat(handles.pitchSyl(1:(expLocal-1)), pitchL, pitchR, handles.pitchSyl((expLocal+1):end));
        handles.fmSyl = horzcat(handles.fmSyl(1:(expLocal-1)), fmL, fmR, handles.fmSyl((expLocal+1):end));
        handles.amSyl = horzcat(handles.amSyl(1:(expLocal-1)), amL, amR, handles.amSyl((expLocal+1):end));
        handles.ampSyl = horzcat(handles.ampSyl(1:(expLocal-1)), ampL, ampR, handles.ampSyl((expLocal+1):end));
        handles.SpecSylinter = cat(3,handles.SpecSylinter(:,:,1:(expLocal-1)), specL, specR, handles.SpecSylinter(:,:,(expLocal+1):end));
        handles.PC1 = horzcat(handles.PC1(1:(expLocal-1)), PC1L, PC1R, handles.PC1((expLocal+1):end));
        handles.PC2 = horzcat(handles.PC2(1:(expLocal-1)), PC2L, PC2R, handles.PC2((expLocal+1):end));
        handles.PC3 = horzcat(handles.PC3(1:(expLocal-1)), PC3L, PC3R, handles.PC3((expLocal+1):end));
    end
end

%Reset ROI stuff
handles = resetROIs(handles);

%Generate indices
[handles.syltempBin, handles.syltempIndx] = createIndices(handles);

handles.syltempStart=handles.sylStart(handles.syltempBin);
handles.syltempEnd=handles.sylEnd(handles.syltempBin);
handles.syltempType=handles.sylType(handles.syltempBin);
handles.syltempIndxKey=handles.IndxKey(handles.syltempBin);

%Redisplay the sorting image
handles=refresh_hist(handles);

%Update the Ungrouped mean syl plot
MeanSpec=log10(mean(handles.SpecSylinter(:,:,:),3));
axes(handles.axes_syl1)
imagesc(MeanSpec)
axis xy

guidata(hObject, handles);

function push_split4_Callback(hObject, eventdata, handles)
%Generate indices
[handles.syltempBin, handles.syltempIndx] = createIndices(handles);

handles.syltempStart=handles.sylStart(handles.syltempBin);
handles.syltempEnd=handles.sylEnd(handles.syltempBin);
handles.syltempType=handles.sylType(handles.syltempBin);
handles.syltempIndxKey=handles.IndxKey(handles.syltempBin);

%Update only those sylTypes that have been defined by the ROIs
%If there is a mask for ROI1, remap sylTypes as specified
if isfield(handles, 'mask3')
    % Use the mask to generate an index of syllable that fall within the ROI
    ROI3Ind = getROIsubset(handles, handles.mask3);
    localPntr = find(ROI3Ind);
    
    % Calculate user-specified splitting threshold in %
    value=get(handles.slider_split4,'Value')/size(handles.MeanSpec(:,:,4),2);

    % Calculate the duration of ROI syllables in seconds based on annotion times
    duration=(handles.syltempEnd(ROI3Ind) - handles.syltempStart(ROI3Ind))/44150;
    margin=100;
    
    % Do the splitting
    %Loop through each syllable in the ROI
    ind = 0;
    a = handles.syltempIndx(ROI3Ind); %creates global pointer array
    for i = 1:length(a)
        j=a(i)+i-1; %Requires sequencing... do not parallelize
        %Calculate the time (in samples) from syllable onset to end of left syllable
        valueaddstart=floor(duration(i)*44150*value+margin);
        
        %Insert new points into the global arrays
        handles.sylStart=horzcat(handles.sylStart(1:j), handles.sylStart(j)+valueaddstart, handles.sylStart(j+1:end));
        handles.sylEnd=horzcat(handles.sylEnd(1:j-1), handles.sylStart(j+1)-floor(2*margin), handles.sylEnd(j:end));
        handles.sylType=vertcat(handles.sylType(1:j-1), str2double(get(handles.edit_split_left4,'String')), str2double(get(handles.edit_split_right4,'String')), handles.sylType(j+1:end));
        handles.IndxKey=vertcat(handles.IndxKey(1:j), handles.IndxKey(j), handles.IndxKey(j+1:end));

        %Update the feature arrays
%         expLocal = localPntr(i) + i - 1;
        expLocal = localPntr(i) + ind;
        ind = ind + 1;
        tempTot = handles.FeatSyl{expLocal};
        snipLength = length(tempTot.Pitch);
        tempL.AM = tempTot.AM(1:floor(snipLength*value)); tempR.AM =  tempTot.AM((floor(snipLength*value)+1):end);
        tempL.FM = tempTot.FM(1:floor(snipLength*value)); tempR.FM =  tempTot.FM((floor(snipLength*value)+1):end);
        tempL.Entropy = tempTot.Entropy(1:floor(snipLength*value)); tempR.Entropy =  tempTot.Entropy((floor(snipLength*value)+1):end);
        tempL.amplitude = tempTot.amplitude(1:floor(snipLength*value)); tempR.amplitude =  tempTot.amplitude((floor(snipLength*value)+1):end);
        tempL.gravity_center = tempTot.gravity_center(1:floor(snipLength*value)); tempR.gravity_center =  tempTot.gravity_center((floor(snipLength*value)+1):end);
        tempL.PitchGoodness = tempTot.PitchGoodness(1:floor(snipLength*value)); tempR.PitchGoodness =  tempTot.PitchGoodness((floor(snipLength*value)+1):end);
        tempL.Pitch = tempTot.Pitch(1:floor(snipLength*value)); tempR.Pitch =  tempTot.Pitch((floor(snipLength*value)+1):end);

        entL = mean(tempL.Entropy); entR = mean(tempR.Entropy);
        pitchL = mean(tempL.Pitch); pitchR = mean(tempR.Pitch);
        fmL = mean(tempL.FM); fmR = mean(tempR.FM);
        amL = mean(tempL.AM); amR = mean(tempR.AM);
        ampL = mean(tempL.amplitude); ampR = mean(tempR.amplitude);
        durL = length(tempL.Entropy)/1000; durR = length(tempR.Entropy)/1000;
        
        %Create PCA projections for the new elements
        fL = [entL, pitchL, fmL, amL, ampL, durL]; fR = [entR, pitchR, fmR, amR, ampR, durR];
        fL_z = (fL-handles.mu)./handles.sigma; fR_z = (fR-handles.mu)./handles.sigma;
        PC1L = fL_z/handles.pcaCoeff(:,1)'; PC1R = fR_z/handles.pcaCoeff(:,1)';
        PC2L = fL_z/handles.pcaCoeff(:,2)'; PC2R = fR_z/handles.pcaCoeff(:,2)';
        PC3L = fL_z/handles.pcaCoeff(:,3)'; PC3R = fR_z/handles.pcaCoeff(:,3)';
        
        %Split the corresponding Spectrogram
        tempSpec = handles.SpecSylinter(:,:,expLocal);
        specDur = size(tempSpec,2);
        specLorig = tempSpec(:,1:(floor(specDur*value)+3)); specRorig = tempSpec(:,(floor(specDur*value)-3):end);
        for k = 1:size(tempSpec,1)
            specL(k,:) = interp1(1:size(specLorig,2), specLorig(k,:), linspace(1, size(specLorig,2), specDur));
            specR(k,:) = interp1(1:size(specRorig,2), specRorig(k,:), linspace(1, size(specRorig,2), specDur));
        end

        if str2double(get(handles.edit_split_left4,'String')) ~= handles.sylvalue
            tempL = []; entL = []; pitchL = []; fmL = []; amL = []; ampL = []; specL = []; PC1L = []; PC2L = []; PC3L = [];
            ind = ind - 1;
        end

        if str2double(get(handles.edit_split_right4,'String')) ~= handles.sylvalue
            tempR = []; entR = []; pitchR = []; fmR = []; amR = []; ampR = []; specR = []; PC1R = []; PC2R = []; PC3R = [];
            ind = ind - 1;
        end

        %Do insertions
        handles.FeatSyl = vertcat(handles.FeatSyl(1:(expLocal-1)), tempL, tempR, handles.FeatSyl((expLocal+1):end));
        handles.entSyl = horzcat(handles.entSyl(1:(expLocal-1)), entL, entR, handles.entSyl((expLocal+1):end));
        handles.pitchSyl = horzcat(handles.pitchSyl(1:(expLocal-1)), pitchL, pitchR, handles.pitchSyl((expLocal+1):end));
        handles.fmSyl = horzcat(handles.fmSyl(1:(expLocal-1)), fmL, fmR, handles.fmSyl((expLocal+1):end));
        handles.amSyl = horzcat(handles.amSyl(1:(expLocal-1)), amL, amR, handles.amSyl((expLocal+1):end));
        handles.ampSyl = horzcat(handles.ampSyl(1:(expLocal-1)), ampL, ampR, handles.ampSyl((expLocal+1):end));
        handles.SpecSylinter = cat(3,handles.SpecSylinter(:,:,1:(expLocal-1)), specL, specR, handles.SpecSylinter(:,:,(expLocal+1):end));
        handles.PC1 = horzcat(handles.PC1(1:(expLocal-1)), PC1L, PC1R, handles.PC1((expLocal+1):end));
        handles.PC2 = horzcat(handles.PC2(1:(expLocal-1)), PC2L, PC2R, handles.PC2((expLocal+1):end));
        handles.PC3 = horzcat(handles.PC3(1:(expLocal-1)), PC3L, PC3R, handles.PC3((expLocal+1):end));
    end
end

%Reset ROI stuff
handles = resetROIs(handles);

%Generate indices
[handles.syltempBin, handles.syltempIndx] = createIndices(handles);

handles.syltempStart=handles.sylStart(handles.syltempBin);
handles.syltempEnd=handles.sylEnd(handles.syltempBin);
handles.syltempType=handles.sylType(handles.syltempBin);
handles.syltempIndxKey=handles.IndxKey(handles.syltempBin);

%Redisplay the sorting image
handles=refresh_hist(handles);

%Update the Ungrouped mean syl plot
MeanSpec=log10(mean(handles.SpecSylinter(:,:,:),3));
axes(handles.axes_syl1)
imagesc(MeanSpec)
axis xy

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bound Section Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_rename1_Callback(hObject, eventdata, handles)

function edit_rename1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_rename2_Callback(hObject, eventdata, handles)

function edit_rename2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_rename3_Callback(hObject, eventdata, handles)

function edit_rename3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_rename4_Callback(hObject, eventdata, handles)

function edit_rename4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function slider_split1_Callback(hObject, eventdata, handles)
axes(handles.axes_syl1)
hold on
imagesc(handles.MeanSpec(:,:,1))
axis xy
plot([get(handles.slider_split1,'Value') get(handles.slider_split1,'Value')], [0.1 100], 'red')
hold off

guidata(hObject, handles);

function slider_split1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function edit_split_left1_Callback(hObject, eventdata, handles)

function edit_split_left1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_split_right1_Callback(hObject, eventdata, handles)

function edit_split_right1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function slider_split2_Callback(hObject, eventdata, handles)
axes(handles.axes_syl2)
hold on
imagesc(handles.MeanSpec(:,:,2))
axis xy
plot([get(handles.slider_split2,'Value') get(handles.slider_split2,'Value')], [0.1 100], 'red')
hold off

guidata(hObject, handles);

function slider_split2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function edit_split_left2_Callback(hObject, eventdata, handles)

function edit_split_left2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_split_right2_Callback(hObject, eventdata, handles)

function edit_split_right2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function slider_split3_Callback(hObject, eventdata, handles)
axes(handles.axes_syl3)
hold on
imagesc(handles.MeanSpec(:,:,3))
axis xy
plot([get(handles.slider_split3,'Value') get(handles.slider_split3,'Value')], [0.1 100], 'red')
hold off

guidata(hObject, handles);

function slider_split3_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function edit_split_left3_Callback(hObject, eventdata, handles)

function edit_split_left3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_split_right3_Callback(hObject, eventdata, handles)

function edit_split_right3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function slider_split4_Callback(hObject, eventdata, handles)
axes(handles.axes_syl4)
hold on
imagesc(handles.MeanSpec(:,:,4))
axis xy
plot([get(handles.slider_split4,'Value') get(handles.slider_split4,'Value')], [0.1 100], 'red')
hold off

guidata(hObject, handles);

function slider_split4_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function edit_split_left4_Callback(hObject, eventdata, handles)

function edit_split_left4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_split_right4_Callback(hObject, eventdata, handles)

function edit_split_right4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
