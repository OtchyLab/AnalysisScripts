function display_motif
%% display motif
Fs=44100;
motif_counter=0;
%% Use system background color for GUI components
panelColor = get(0,'DefaultUicontrolBackgroundColor');
global folder_dir;
global time;
global hl;
global last;
global htext;
counter=1;
[filename,dirpath]=uigetfile('*.*');
old_dir=cd;
cd (dirpath);
all_files=dir(dirpath);
[file_order,file_no]=sortfiles(all_files, filename);
n=file_no;
filename=all_files(file_order(n-2)).name

%% ------------ GUI layout ---------------

%% Set up the figure and defaults
f = figure('Units','characters',...
        'Position',[2 30 200 35],...
        'Color',panelColor,...
        'IntegerHandle','off',...
        'Renderer','painters',...
        'Toolbar','figure',...
        'NumberTitle','off',...
        'Name','Spectrogram Display',...
        'ResizeFcn',@figResize,...
        'CloseRequestFcn',@figClose,...
        'KeyPressFcn',@keypress);

%% Create the bottom uipanel
botPanel = uipanel('BackgroundColor',panelColor,...
    'Units','characters',...
    'Position',[1/20 1/20 200 8],...
    'Parent',f,...
    'ResizeFcn',@botPanelResize);


%% Create the center panel
centerPanel = uipanel('bordertype','etchedin',...
    'Units','characters',...
    'Position', [1/20 8 180 27],...
    'Parent',f);

%% Create the right panel
rightPanel = uipanel('bordertype','etchedin',...
    'Units','characters',...
    'Position', [180 8 20 27],...
    'Parent',f,...
    'ResizeFcn',@rightPanelResize);

%% Add an axes to the center panel
a = axes('parent',centerPanel);
%% Add buttons      
playButton = uicontrol(f,'Style','pushbutton','Units','characters',...
        'Position',[45 1 24 2],...
        'String','Play Song',...
        'Parent',botPanel,...
        'Callback',@songplay);
previousButton = uicontrol(f,'Style','pushbutton','Units','characters',...
        'Position',[10 1 24 2],...
        'String','previous Song',...
        'Parent',botPanel,...
        'Callback',@previous_song);
nextButton = uicontrol(f,'Style','pushbutton','Units','characters',...
        'Position',[80 1 24 2],...
        'String','Next Song',...
        'BackgroundColor',panelColor,...
        'Parent',botPanel,...
        'Callback',@next_song);
fileSlider = uicontrol(f,'Style','slider','Units','characters','Min',1,'Max',length(all_files)-2,...
        'Position',[20 4 75 1],...
        'BackgroundColor',panelColor,...
        'String',n-2,...
        'SliderStep',[1/(length(all_files)-3) 2/(length(all_files)-3)],...
        'Value',n-2,...
        'Parent',botPanel,...
        'Callback',@sliderpos);
fileBox = uicontrol(f,'Style','text','Units','characters',...
        'Position',[20 5.5 75 2],...
        'BackgroundColor',panelColor,...
        'String',filename,...
        'Parent',botPanel);
    
folderButton = uicontrol(f,'Style','pushbutton','Units','characters',...
        'Position',[10 7 18 2],...
        'String','Save',...
        'BackgroundColor',panelColor,...
        'Parent',rightPanel,...
        'Callback',@new_folder);


counterButton = uicontrol(f,'Style','text','Units','characters',...
        'Position',[25 6 3 1],...
        'String',num2str(motif_counter),...
        'BackgroundColor',panelColor,...
        'Parent',rightPanel);
   % 'Callback',{@syllfnc,i}
motifButton = uicontrol(f,'Style','pushbutton','Units','characters',...
        'Position',[10 6 10 2],...
        'String','Motif',...
        'BackgroundColor',panelColor,...
        'Parent',rightPanel,...
        'Callback',{@motiffnc});

%display the spectrogram of current song


%% ------------ Callback Functions ---------------


function figClose(src,evt)
    cd(old_dir);
    closereq;
    motif_counter=0;
end
% Figure resize function
function figResize(src,evt)
    fpos = get(f,'Position');
    set(botPanel,'Position',...
        [1/20 1/20 fpos(3)-.1 fpos(4)*8/35])
    set(centerPanel,'Position',...
        [1/20 fpos(4)*8/35 fpos(3)*180/200 fpos(4)*27/35]);
    set(rightPanel,'Position',...
        [fpos(3)*180/200 fpos(4)*8/35 fpos(3)*20/200 fpos(4)*27/35]);
end
    
% Bottom panel resize function
function botPanelResize(src, evt)
    bpos = get(botPanel,'Position');
    set(previousButton,'Position',...
        [bpos(3)*10/120 bpos(4)*1/8 bpos(3)*24/120 2])
    set(playButton,'Position',...
        [bpos(3)*45/120 bpos(4)*1/8 bpos(3)*24/120 2])
    set(nextButton,'Position',...
        [bpos(3)*80/120 bpos(4)*1/8 bpos(3)*24/120 2])
    set(fileSlider,'Position',...
        [bpos(3)*20/120 bpos(4)*4/8 bpos(3)*75/120 1])
    set(fileBox,'Position',...
        [bpos(3)*20/120 bpos(4)*5.5/8 bpos(3)*75/120 2])

end

% Right panel resize function
function rightPanelResize(src,evt)
    rpos = get(rightPanel,'Position');
    set(folderButton,'Position',...
        [rpos(3)*2/20 rpos(4)*7.3/8 rpos(3)*16/20 2])
    set(counterButton,'Position',...
        [rpos(3)*14/20 rpos(4)*(6.3)/8 rpos(3)*3/20 1]);
    set(motifButton,'Position',...
        [rpos(3)*2/20 rpos(4)*(6.3)/8 rpos(3)*10/20 2]);
 
end

%% Callback for list box
    
function songplay(obj,event)
filename=all_files(file_order(n-2)).name
song=wavread(filename);
wavplay(song,Fs);
end
function next_song(obj,event)
    
    if (n<length(all_files))
        n=n+1;
        cd(dirpath);
        [song,filename]=disp_song(n, all_files,file_order);
        spectrogram = findobj('Type','image');
        set(spectrogram,'ButtonDownFcn',@mousepress);
        set(fileSlider,'Value',n-2);
        set(fileBox,'String',filename);
        counter=1;
        
        
    else
        
    end
end
function previous_song(obj,event)
    
    if (n>3)
        n=n-1;
        cd(dirpath);
        [song,filename]=disp_song(n, all_files,file_order);
        spectrogram = findobj('Type','image');
        set(spectrogram,'ButtonDownFcn',@mousepress);
        set(fileSlider,'Value',n-2);
        set(fileBox,'String',filename);
        counter=1;
        
    end
end
function sliderpos(obj,event)
    cd(dirpath);
    n=floor(get(fileSlider, 'Value'))+2
    [song,filename]=disp_song(n, all_files, file_order);
    spectrogram = findobj('Type','image');
    set(spectrogram,'ButtonDownFcn',@mousepress);
    set(fileBox,'String',filename);
    counter=1;
end
function motiffnc(obj,event)
    %read in the coordinates of the mouse and save the part of the file in
    %the appropriate folder.
    old_dir=cd;
    new_dir=[folder_dir];
    cd (new_dir);
    motiftosave=song(round(min(time)*Fs):round(max(time)*Fs),:);
    wavwrite(motiftosave,Fs,[num2str(min(time)),'_',num2str(max(time)),filename]);
    cd (old_dir);
    counter=1;
    motif_counter=motif_counter+1;
    htext=text(mean(time),8000,'motif','HorizontalAlignment','center','Color',[1 1 0]);
    set(counterButton,'String',num2str(motif_counter));
    delete(hl);
 
    
end

function new_folder(obj,event)
    %here a dialog box to choose a folder name
    folder = inputdlg({'Folder to save the syllables','How many syllables?'});
    mkdir(folder{1});
    folder_dir=[dirpath,folder{1}]; %folder for the syllables
    cd (folder_dir);
    mkdir('Motifs'); 
    cd(dirpath);
end



function mousepress(obj,event)   
    if (counter<3)    
        click = get(a,'CurrentPoint');
        time(counter)=click(1,1);
        hl(counter) = line([time(counter) time(counter)],[0 10000]);
        set(hl(counter),'Color',[1 1 1]);
    else
        ncounter=abs(mod(counter,2)-2);
        delete(hl(ncounter)); 
        click = get(a,'CurrentPoint');
        time(ncounter)=click(1,1);
        hl(ncounter) = line([time(ncounter) time(ncounter)],[0 10000]);
        set(hl(ncounter),'Color',[1 1 1]);
    end
counter=counter+1;
end

axes(a);
song=wavread(filename);
specgram1(song,512,44100,400,350);ylim([0 10000]);
spectrogram = findobj('Type','image');
set(spectrogram,'ButtonDownFcn',@mousepress);



end % uipanel1





