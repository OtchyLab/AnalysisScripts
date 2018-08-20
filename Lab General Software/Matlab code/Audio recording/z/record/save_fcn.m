userdata.directed=directed;
userdata.dirname=dirname;
userdata.filecount=filecount;
if strcmp(caller,'Feedback')
    userdata.soft_trig_time=soft_trig_time;
    userdata.golay_times=golay_times;
    if do_feedback
        %if noise_starts(1)==0;
        %    noise_starts(1)=soft_trig_time*scanrate;
        %end
    end
    userdata.noise_starts=noise_starts;
    userdata.noise_stops=noise_stops;
    userdata.do_feedback=do_feedback;
    userdata.on_amps=on_amps;
    userdata.off_amps=off_amps;
    userdata.samp_times=samp_times;
    userdata.thr_amp=thr_amp;
end
save(filename_userdata,'userdata');
%%%%%%%%%%%%%%%%% check for new date
if strcmp(dirname,datestr(date,29))
    namehelp=NextName(namehelp);
else
    dirname_new=datestr(date,29);
    set(hDirectory,'string',dirname_new);
    if ~(exist([savepath birdname '/' dirname_new])==7)
        mkdir([savepath birdname],dirname_new);
    end
    eval(['!move ' savepath_win birdname '\' dirname '\' namehelp '.daq ' savepath_win birdname '\' dirname_new '\' birdname songtype_prefix '-' caller_prefix '0001.daq']);
    eval(['!move ' savepath_win birdname '\' dirname '\U_' namehelp '.mat ' savepath_win birdname '\' dirname_new '\U_' birdname songtype_prefix '-' caller_prefix '0001.mat']);
    dirname=dirname_new;
    filecount='0002';
    namehelp=[birdname songtype_prefix '-' caller_prefix filecount]; 
end
set(hFilecount,'string',namehelp(end-3:end));