if state_ai==0 & ~(do_record & ~is_record)
  do_record=~do_record;
set(htoggle_record,'backgroundcolor',[~do_record do_record 0]); % 
  if do_record
  set(ai, 'TimerFcn', '');
set(hPatch,'facecolor',[.3 0.3 1]);drawnow

         birdname=get(hBirdname,'string');
        if ~(exist([savepath birdname])==7)
            mkdir(savepath,birdname);
        end
        if ~(exist([savepath birdname '/' dirname])==7)
            mkdir([savepath birdname],dirname);
        end
     
        caller_prefix=set_caller_prefix(caller,0); % no feedback
        filecount=get(hFilecount,'string');
        namehelpS=[birdname  'S' songtype_prefix '-'  filecount '-' datestr(now,30)];    % 'S' for spontaneous
        filename=[savepath birdname '/' dirname '/' namehelpS '.daq'];
        filename_userdata=[savepath birdname '/' dirname '/U_' namehelpS '.mat'];
        set(ai,'LoggingMode','Disk&Memory');
        set(ai,'LogToDiskMode','index');
        set(ai,'LogFileName',filename);    
    start(ai);   % off we go !!
pause(pre_trigger_int);
trigger(ai);
is_record=1;
else
    y=getdata(ai,ai.SamplesAvailable,'native')';
    set(hPatch,'facecolor',[0 0 0]);
    for i=1:length(ai.channel)
        maxi=max(abs(double(y(i,:))));
        if maxi==2^nbits
            warndlg(['Clipping on ' ai.Channel.ChannelName(i) ', in file ' filename]);
            end
        set(hPatch_el(i),'xdata',ind_el_height/nbits/log(2)*log(maxi)*[0 0 1 1]);
    end
        stop(ai);flushdata(ai);   
  userdata.directed=directed;
userdata.dirname=dirname;
userdata.filecount=filecount;
  save(filename_userdata,'userdata');
    stop(ai); is_record=0;
end
end