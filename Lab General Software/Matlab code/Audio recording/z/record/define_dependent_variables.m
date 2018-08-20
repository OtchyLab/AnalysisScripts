%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DEFINE SOME DEPENDENT VARIABLES
path(path,'..');path(path,'../analyze');
state_ai=0;  % 1=running    0=stop
loudspeaker_warning=0;

if ~exist('do_save')
    do_save=do_save_default;
end
if ~exist('do_plot')
    do_plot=do_plot_default;
end
if ~exist('directed')
    directed=directed_default;
end
do_record=0; is_record=0;
i=regexp(savepath,'\/'); savepath_win=savepath; savepath_win(i)='\';
dirname= datestr(date,29);  % directory for saving of data
pre_trigger_samples=ceil(pre_trigger_int*scanrate);
mid_trigger_samples=ceil(mid_trigger_int*scanrate);
post_trigger_samples=ceil(post_trigger_int*scanrate);
peek_samples=peek_period*scanrate; %
%trigger_cutout_i=round(on_mic/scanrate*trigger_samples);

if strcmp(caller,'Feedback');

%%%%%%% NOISE
noise_samples=noise_int*scanrate_o; % duration of noise pattern in samples, a power of 2
noise_cutout_i=round(cutout_noise/scanrate_o*noise_samples);

noise_envelope=ones(noise_samples,1);
ramp_samples=floor(noise_samples/10);
%noise_envelope(1:ramp_samples)=1./(1+exp(-(((1:ramp_samples)-ramp_samples/2)/ramp_samples*12)));
noise_envelope(1:ramp_samples)=(0:ramp_samples-1)/(ramp_samples-1);
noise_envelope(end:-1:end-ramp_samples+1)=noise_envelope(1:ramp_samples);

    
    if ~exist('pfeedback')
    pfeedback=pfeedback_default;
end
    if on_mic(2)>cutout_noise(2) | on_mic(1)<cutout_noise(1)
    fprintf('Problem with noise-microphone cutout \n');
    break
end
if pre_golay_int>mid_trigger_int | pre_golay_int<noise_int
    fprintf('Problem with onset of Golay pulses\n');
    break
end
end

thr_amp=3;% trigger_threshold (scale), dot not change this value
max_amp=2*thr_amp; % maximum amplitude for display
axis_amp=2*thr_amp; % axis limits for display
on_amp=0;  % current ON-amplitude
off_amp=0; % current OFF-amplitude

bytes=8000;
userdata.caller=caller;
userdata.birdname=birdname;
userdata.micro=micro;
userdata.speaker_input=speaker_input;
userdata.channels=channels;
%userdata.chans=chans;
userdata.pre_trigger_samples=pre_trigger_samples;
userdata.mid_trigger_samples=mid_trigger_samples;
userdata.post_trigger_samples=post_trigger_samples;

if strcmp(caller,'Feedback');
    userdata.channels_o=channels_o;
    userdata.scanrate_o=scanrate_o;
    userdata.speaker_input=speaker_input;
    userdata.cutout_noise=cutout_noise;
    userdata.on_mic=on_mic;
    userdata.intergolay_int=intergolay_int;
end
 