%%%%%%%%%%% NEW NOISE
w_noise=1-2*rand(noise_samples,1); 
%a=fft(w_noise);a(noise_cutout_i(1):noise_cutout_i(2))=0; 
%w_noise=abs(ifft(a));
w_noise=noise_envelope.*filter(filt_noise_b,filt_noise_a,w_noise);
noise_starts=[]; % times of noise onsets
noise_stops=[]; % times of noise offsets

%%%%%%%%%%% NEW GOLAY
[golay_dat1 golay_dat2]=golay(golay_order);
golay_dat=2*[golay_dat1 zeros(1,intergolay_int*scanrate_o) golay_dat2 0]';
length_golay_dat=length(golay_dat);
golay_up=0;  % are we ready to play a golay-pulse?
golay_times=[]; % times of golay pulse onsets 

userdata.w_noise=w_noise;
userdata.golay_dat=golay_dat;