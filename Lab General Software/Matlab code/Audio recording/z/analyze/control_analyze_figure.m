%%%%%%%%%%%%%%%%%% Get synopsis
if exist([loadpath birdnameA '/synopsis.mat'])==2
    load([loadpath birdnameA '/synopsis.mat'])
else
    synopsis='No synopsis';
end

% Position the GUI in the middle of the screen
btnColor=get(0,'DefaultUIControlBackgroundColor');
screenUnits=get(0,'Units');
set(0,'Units','pixels');
screenSize=get(0,'ScreenSize');
set(0,'Units',screenUnits); 
figWidthA=screenSize(3); figHeightA=screenSize(4)/2; %size of main screen
figPosA=[(screenSize(3)-figWidthA)/2 30         figWidthA                    figHeightA+30];%position of the main figure on the main screen
ctr_sizeA=300;
figPos_ctrA=[screenSize(3)-2*ctr_sizeA-8 screenSize(4)-ctr_sizeA ctr_sizeA ctr_sizeA-20]; %position for the control window
% Create the figure window.

% Create the control figure 
hFig_ctrA=figure(...                    
    'Color'             ,btnColor                 ,...
    'IntegerHandle'     ,'on'                    ,...
    'DoubleBuffer'      ,'on'                     ,...
    'DeleteFcn'         ,'close_analyze_window',...
    'MenuBar'           ,'none'                   ,...
    'HandleVisibility'  ,'on'                     ,...
    'Name'              ,'Control: Analyze'  ,...
    'Tag'               ,'ControlA'  ,...
    'NumberTitle'       ,'off'                    ,...
    'Units'             ,'pixels'                 ,...
    'Position'          ,figPos_ctrA                   ,...
    'UserData'          ,[]                       ,...
    'Colormap'          ,[]                       ,...
    'Pointer'           ,'arrow'                  ,...
    'Visible'           ,'off'                     ...
    );


hAxes_ctrA = axes(... 
    'Position'          , [0.0200 0.0200 0.98 0.98],...
    'Parent'            , hFig_ctrA,...
    'XLim'              , [0 1]...
    );
hText_ctrA = text(...
            'Parent',hAxes_ctrA,...
            'String','','position',[0.02, 0.1 ],...
            'FontSize',14);
%%%%%%%%%%%%%%% Create the pushbuttons
htoggle_plotA = uicontrol(...
    'Parent'          , hFig_ctrA,...
    'Style'           , 'pushbutton',...
    'Units'           , 'normalized',...
    'Position'        , [0.02 0.86 0.2 0.12],...
    'Value'           , 1,...
    'backgroundcolor' , [~do_plotA do_plotA 0],...%switches between red and green
    'String'          , 'PLOT',...
    'fontsize'        , 14,...
    'Callback'        , 'do_plotA=~do_plotA;set(htoggle_plotA,''backgroundcolor'',[~do_plotA do_plotA 0]); if do_plotA; control_analyze_figure_plot;analyze_plot;else; delete(hFigA); end;');

%%%%%%%%%%%%%%%%%%% create the other buttons etc
hBirdnameA=uicontrol('Parent',hFig_ctrA,'Style','popupmenu','String',birdnames,...
    'units','normalized','position',[0.02 0.74 0.4 0.1],'FontSize',14,...
    'HorizontalAlignment','left','Value',length(birdnames),...
    'Callback','birdnameA=birdnames{get(hBirdnameA,''Value'')};[dirnames,dirnameA]=get_directories(loadpath,birdnameA);set(hDirectoryA,''string'',dirnames);get_filenames;analyze_plot');

hDirectoryA=uicontrol('Parent',hFig_ctrA,'Style','popupmenu','String',dirnames,...
    'units','normalized','position',[0.55 0.74 0.40 0.1],'FontSize',14,...
    'HorizontalAlignment','left','Value',length(dirnames),...
    'Callback','dirnameA=dirnames{get(hDirectoryA,''Value'')};get_filenames;analyze_plot');

hSongtypeA=uicontrol('Parent',hFig_ctrA,'Style','popupmenu','String',callers.songtype_text{caller_iA},...
    'units','normalized','position',[0.02 0.6 0.45 0.065],'fontsize',14,...
    'HorizontalAlignment','left','Value',get_index(callers.songtype_text{caller_iA},callerA.songtype_text),...
    'Callback','callerA.songtype=callers.songtype{caller_iA}(get(hSongtypeA,''Value''));callerA.songtype_text=callers.songtype_text{caller_iA}(get(hSongtypeA,''Value''));get_filenames;analyze_plot');

hCaller=uicontrol('Parent',hFig_ctrA,'Style','popupmenu','String',callers.prefix_text{caller_iA},...
    'units','normalized','position',[0.55 0.6 0.45 0.065],'FontSize',14,...
    'HorizontalAlignment','left','Value',get_index(callers.prefix_text{caller_iA},callerA.prefix_text),...
    'Callback','callerA.prefix=callers.prefix{caller_iA}(get(hCaller,''Value''));callerA.prefix_text=callers.prefix_text{caller_iA}(get(hCaller,''Value''));get_filenames;analyze_plot');

hFilecountA=uicontrol('Parent',hFig_ctrA,'Style','popupmenu','String',filecounts,...
    'units','normalized','position',[0.55 0.44 0.24 0.065],'FontSize',14,...
    'HorizontalAlignment','left','Value',filecountA,...
    'Callback','filecountA=get(hFilecountA,''Value'');namehelpA=namehelpsA{filecountA};analyze_plot');

hFilecountUpA=uicontrol('Parent',hFig_ctrA,'Style','pushbutton','String','Prev',...
    'units','normalized','position',[0.81 0.45 0.15 0.065],'FontSize',14,... 
    'Callback','filecountA=filecountA-1;namehelpA=namehelpsA{filecountA};analyze_plot');

hFilecountDownA=uicontrol('Parent',hFig_ctrA,'Style','pushbutton','String','Next',...
    'units','normalized','position',[0.81 0.37 0.15 0.065],'FontSize',14,... 
    'Callback','if filecountA<length(filecounts);filecountA=filecountA+1;namehelpA=namehelpsA{filecountA};analyze_plot;end;');

hPlaySong=uicontrol('Parent',hFig_ctrA,'Style','pushbutton','String','Play',...
    'units','normalized','position',[0.02 0.40 0.20 0.1],'FontSize',14,... 
    'Callback','soundsc(yA,scanrateA)');

hWavefile=uicontrol('Parent',hFig_ctrA,'Style','pushbutton','String','Wav',...
    'units','normalized','position',[0.32 0.40 0.20 0.1],'FontSize',14,... 
    'Callback','wave_save');


hOnline=uicontrol('Parent',hFig_ctrA,'Style','radiobutton','String','Online',...
    'units','normalized','position',[0.36 0.9 0.3 0.065],'fontsize',14,...
    'buttondownfcn','',...
    'callback','if ~online; online=1;online_vs_offline;end;',...
    'Value',online);
hOffline=uicontrol('Parent',hFig_ctrA,'Style','radiobutton','String',...
    'Offline','units','normalized','position',[0.7 0.9 0.4 0.065],...
    'fontsize',14,'callback','if online online=0;online_vs_offline;end;',...
    'buttondownfcn','',...
    'Value',~online,'HandleVisibility','off');
hSynopsisA=uicontrol('Parent',hFig_ctrA,'Style','text','String',synopsis,...
    'units','normalized','position',[0.02 0.24 0.98 0.12],'FontSize',11,...
    'HorizontalAlignment','left','backgroundcolor',0*[1 1 1],'foregroundcolor',1*[1 1 1]);



%  hmenu(1) = uimenu('Parent', hFig,'Label', 'File');
%  hmenu(2) = uimenu(hmenu(1),...
%      'Label', 'Close demoai_fft',...
%      'Callback', 'demoai_fft(''close'',gcbf)');
%  hmenu(3) = uimenu('Parent', hFig,...
%      'Label', 'Help');
%  hmenu(4) = uimenu(hmenu(3),...
%      'Label', 'Data Acquisition Toolbox',...
%      'Callback', 'helpwin(''daq'')');
%  hmenu(5) = uimenu(hmenu(3),...
%      'Label', 'demoai_fft',...
%      'Callback', 'helpwin(''demoai_fft'')');
axis off;

%    handle_ctr.filename=hfilename;

set(hAxes_ctrA, 'HandleVisibility', 'on');
set(hFig_ctrA,'Visible','on','HandleVisibility', 'on');
zoom off
