function [warped] = warpit(d1,d2,sr,p,q)
%Dynammic Time Warping for one time series, d2, with parameters previouly
%found with dtw.

 % Constants for FFT
 nfft = 256;
 window = 256;
 overlap = window * .75;

 % Calculate STFT features for sound (25% window overlap)
 D1 = spectrogram(d1,window,overlap,nfft,sr);
 D2 = spectrogram(d2,window,overlap,nfft,sr);

 % Calculate the frames in D2 that are indicated to match each frame
 % in D1, so we can resynthesize a warped, aligned version
 D2i1 = zeros(1, size(D1,2));
 for i = 1:length(D2i1); 
    D2i1(i) = q(min(find(p >= i))); 
 end
 
 % Phase-vocoder interpolate D2's STFT under the time warp
 D2x = pvsample(D2, D2i1-1, nfft/2);
 
 % Invert it back to time domain
 %d2x = istft(D2x, 512, 512, 128);
 d2x = istft(D2x, nfft, window, window-overlap);

 % Warped version added to original target (have to fine-tune length)
 warped = resize(d2x', length(d1),1);
 
end
