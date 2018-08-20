daqreset;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ANALOG INPUT
ai = analoginput(input_device,input_hwaddress);
if strcmp(input_device,'nidaq')
    set(ai,'InputType','SingleEnded');
    %set(ai,'TransferMode','Interrupts');
end

channels=sort(channels);
chans_col={};  % define channel list and colors
if micro>-1
    addchannel(ai,micro,'Microphone'); chans_col={'k'};
end 
if speaker_input>-1
    addchannel(ai,speaker_input,'Loudspeaker');
    chans_col{end+1}='r';
end
if length(channels)>0
    for i=1:length(channels)
        addchannel(ai,channels(i),['Channel ' num2str(channels(i))]);
        chans_col{end+1}='b';
    end
end
set(ai,'RuntimeErrorFcn','disp(''There is a Problem. You may want to consider increasing the buffer size for the analog input.'')');
set(ai,'StartFcn','t_start=clock;');
set(ai, 'TriggerType', 'manual');
set(ai, 'TriggerDelayUnits', 'seconds');
set(ai, 'SampleRate', scanrate);
if strcmp(caller,'Feedback')
    set(ai, 'SamplesPerTrigger',inf);
%ti=timer;
%set(ti,'BusyMode','drop');
%set(ti,'Period',rep_fcn_period);
%set(ti,'ExecutionMode','fixedRate');
%set(ti,'TasksToExecute',100);

    set(ai, 'TimerPeriod', rep_fcn_period);  
    set(ai,'ManualTriggerHwOn','start');
    set(ai, 'TriggerRepeat', 0);
    set(ai, 'TriggerDelay', -pre_trigger_int);
elseif strcmp(caller,'Playback')
    set(ai, 'SamplesPerTrigger',bos_samples);
    set(ai,'ManualTriggerHwOn','start');
    set(ai,'TriggerRepeat',1);
end
set(ai,'BufferingConfig',[4096 2000]);
%set(ai,'InputOverRangeFcn','warndlg(''Clipping'')');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ANALOG OUTPUT
ao=analogoutput(output_device,output_hwaddress);
addchannel(ao,channels_o,'Loudspeaker');
set(ao,'SampleRate',scanrate_o);
set(ao,'TriggerType','Manual');
putsample(ao,0);
if strcmp(caller,'Feedback')
    if do_off_trig
            set(ao, 'RepeatOutput', 0);
%                        set(ao, 'RepeatOutput', inf); % for rep_fun_continuous
  %  set(ao,'StopFcn','putsample(ao,0);');
        else
        %set(ao, 'RepeatOutput', inf);
        set(ao, 'RepeatOutput', 0);
end
elseif strcmp(caller,'Playback')
    set(ao, 'RepeatOutput', 0);
end    
