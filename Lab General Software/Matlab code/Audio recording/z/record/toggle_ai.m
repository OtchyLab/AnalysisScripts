%if strcmp(ai.Running, 'On') 
if state_ai==1
    % Stop the device.
    if ai.TriggersExecuted==0  
        set(ai, 'TimerFcn', '');
        %stop(ti); 
        stop(ai); stop(ao);
    end
    %set(handle_ctr.toggle_start, 'String', 'RUN');
    set(htoggle_start,'backgroundcolor',[1 0 0]); % red
    % Store the new state.
    state_ai = 0;
elseif do_record==0
    % Toggle the Start/Stop string.
    %set(handle_ctr.toggle_start, 'String', 'STOP');
    set(htoggle_start,'backgroundcolor',[0 1 0]) % green
    % Store the new state.
    state_ai = 1;
        %if strcmp(caller,'Feedback')
    pfeedback=max(0,min(1,str2num(get(hPfeedback,'string'))));        
    if isempty(pfeedback)
        pfeedback=pfeedback_default;
    end
    if pfeedback>0 & speaker_input==-1 & loudspeaker_warning==0
        loudspeaker_warning=1;
        warndlg('You should record the Loudspeaker output');
    end
    set(hPfeedback,'string',num2str(pfeedback));
    do_feedback=rand(1)<pfeedback;
    
    
    %%%%%%%%%% PRE-SAVE 
    if do_save  
        birdname=get(hBirdname,'string');
        if ~(exist([savepath birdname])==7)
            mkdir(savepath,birdname);
        end
        if ~(exist([savepath birdname '/' dirname])==7)
            mkdir([savepath birdname],dirname);
        end
        songtype_prefix=set_songtype_prefix(caller,get(hDirected,'Value'));
        caller_prefix=set_caller_prefix(caller,do_feedback);
        filecount=get(hFilecount,'string');
        namehelp=[birdname  songtype_prefix '-' caller_prefix filecount];
        filename=[savepath birdname '/' dirname '/' namehelp '.daq'];
        filename_userdata=[savepath birdname '/' dirname '/U_' namehelp '.mat'];
        set(ai,'LoggingMode','Disk&Memory');
        set(ai,'LogToDiskMode','index');
        set(ai,'LogFileName',filename);    
    else
        set(ai,'LoggingMode','Memory');
        set(ai,'LogFileName','');
    end
    
    % Start the device.    
    %tic;
    %flushdata(ai);
    %toc;
    start(ai); still_trig_low=clock;
    pause(peek_period+0.1);
    soft_trig=1; real_trig=0; post_trig=0; stop_trig=0;
   putdata(ao,w_noise);
   start(ao);
   noise_up=0; samp_times=[];on_amps=[];off_amps=[];
   ai.TimerFcn='rep_fun'; 
end
