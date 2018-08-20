%t_now=etime(clock,t_now)
%t_now=clock;
%if state_ai==1
%%%%%%%%%%%%%%%% PEEKDATA     
x = peekdata(ai, peek_samples);
x_on=filter(on_mic_b,on_mic_a,x(:,1));
on_amp=log(1+on_amp_scale*norm(x_on)^2);
if do_off_trig
    x_off=filter(off_mic_b,off_mic_a,x(:,1));
    off_amp=log(1+off_amp_scale*norm(x_off)^2);
end
set(hLine_ctr,'xdata',on_amp*[1 1]); %drawnow

if  on_amp>thr_amp & do_off_trig*off_amp<thr_amp
    still_trig_low=clock;
    t_still=0;
else
    t_still=etime(clock,still_trig_low);
end


%%%% SOFT TRIG
if  soft_trig & on_amp>thr_amp & do_off_trig*off_amp<thr_amp   % a trigger event occured
    soft_trig=0; real_trig=1;
    soft_trig_time_abs=clock;
    soft_trig_time=min(pre_trigger_int,etime(clock,t_start));
    set(hPatch,'facecolor',[.3 0.3 1]);drawnow
end


%if ai.TriggersExecuted>0
 if ~soft_trig
     samp_times=[samp_times soft_trig_time+etime(clock,soft_trig_time_abs)];
    on_amps=[on_amps on_amp];
    off_amps=[off_amps off_amp];
end

%%%%%%%%%%%%%%%%%% NOISE and GOLAY
if do_feedback 
    if  strcmp(ao.Sending,'Off')
        if on_amp>thr_amp & do_off_trig*off_amp<thr_amp %%% A trigger event
            if ~golay_up | (golay_up  & t_still<pre_golay_int) % still noise
              if ao.SamplesAvailable==0
                  putdata(ao,w_noise);
              end
                  if strcmp(ao.Running,'Off')
                   start(ao);
               end
               trigger(ao);
             noise_up=1;
               if ai.SamplesAvailable==0
                    noise_starts=[noise_starts (soft_trig_time+etime(clock,soft_trig_time_abs))*ai.SampleRate];
                else
                    noise_starts=[noise_starts get(ai,'SamplesAvailable')];
                end
                golay_up=1;    
            end       
            % no trigger-event  
        elseif   golay_up & t_still>pre_golay_int     % golay  
            stop(ao);
          %  set(ao,'RepeatOutput',0);
            putdata(ao,golay_dat); 
            start(ao); 
            trigger(ao);
            golay_times=[golay_times get(ai,'SamplesAvailable')];
            ai.TimerFcn='';
            pause(length_golay_dat/scanrate_o+wait_post_golay);
            golay_up=0; stop(ao); %set(ao,'RepeatOutput',inf);            
            %putdata(ao,w_noise);
            ai.TimerFcn='rep_fun'; 
        end
    end
    if  strcmp(ao.Running,'Off') & noise_up & t_still>noise_int  % Running= 'On', stop noise
        % while rem(ao.SamplesOutput+ramp_samples/4,noise_samples)>ramp_samples/2; end
        ai.TimerFcn='';
        % stop(ao);
        noise_up=0;
        putsample(ao,0);putdata(ao,w_noise);
        start(ao);
        ai.TimerFcn='rep_fun';
        if ai.SamplesAvailable==0
            noise_stops=[noise_stops (soft_trig_time+etime(clock,soft_trig_time_abs))*ai.SampleRate];
        else
            noise_stops=[noise_stops get(ai,'SamplesAvailable')];
        end
    end
end

%%%%% REAL TRIG    
if real_trig & etime(clock,t_start)>pre_trigger_int+.01
    trigger(ai);
    real_trig=0; post_trig=1;  
    %if ai.TriggersExecuted==1 & etime(clock,still_trig_low)>pre_golay_int & golay_up
    
    %%%%% POST TRIG    
elseif post_trig & t_still>mid_trigger_int
    post_trig=0; stop_trig=1;
    stop_onset=clock;       
    if do_plot
        pre_plot_time=clock;
        ai.TimerFcn='';
        plot_before_post; drawnow;
        if etime(clock,pre_plot_time)>post_trigger_int
            warndlg(['Plotting takes ' num2str(etime(clock,pre_plot_time)/post_trigger_int,3) ' times as long as the post-trigger interval']);
        end
        ai.TimerFcn='rep_fun';
    end    
    %%%%% STOP TRIG    
elseif stop_trig & etime(clock,stop_onset)>post_trigger_int
    stop_trig=0;
    ai.TimerFcn='';     
    y=getdata(ai,ai.SamplesAvailable,'native')';
    set(hPatch,'facecolor',[0 0 0]);
    for i=1:length(ai.channel)
        maxi=max(abs(double(y(i,:))));
        if maxi==2^nbits
            if do_save
            warndlg(['Clipping on ' ai.Channel.ChannelName(i) ', in file ' filename]);
        else
            warndlg(['Clipping on ' ai.Channel.ChannelName(i)]);
        end
        
        end
        set(hPatch_el(i),'xdata',ind_el_height/nbits/log(2)*log(maxi)*[0 0 1 1]);
    end
        stop(ai);flushdata(ai);   stop(ao);
    if do_save 
        save_fcn;
    end
    if exist('online') & exist('do_plotA')
        if online & do_plotA
            % while strcmp(ai.Running,'On')
            % end;
            do_feedback_A=do_feedback;
            userdata.noise_starts=noise_starts;
            userdata.noise_stops=noise_stops;
            userdata.golay_times=golay_times;         
            userdata.post_trigger_samples=post_trigger_samples;
            userdata.golay_dat=golay_dat;
            analyze_plot;drawnow;
        end
    end
    if strcmp(caller,'Feedback')
        new_noise_golay;
    end
    
    if state_ai==1
        state_ai=0;
        toggle_ai; 
    end
end
%end
