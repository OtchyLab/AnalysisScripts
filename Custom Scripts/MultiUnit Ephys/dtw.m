function [d2x,p,q] = dtw(d1,d2,sr)
%Dynammic Time Warping for one time series, d2, onto the structure of a
%second time series, d1.  Both time series should be sampled at a single
%frequency, sr -- though this constraint can be removed in the first line
%of code.

 % Constants for FFT
 nfft = 256;
 window = 256;
 overlap = window * .75;

 % Calculate STFT features for both sounds (25% window overlap)
 D1 = spectrogram(d1,window,overlap,nfft,sr);
 D2 = spectrogram(d2,window,overlap,nfft,sr);

 % Construct the 'local match' scores matrix as the cosine distance 
 % between the STFT magnitudes
 SM = simmx(abs(D1),abs(D2));
 
 % Look at it:
%  subplot(121)
%  imagesc(SM)
%  colormap(1-gray)
 
 % You can see a dark stripe (high similarity values) approximately
 % down the leading diagonal.

 % Use dynamic programming to find the lowest-cost path between the 
 % opposite corners of the cost matrix
 % Note that we use 1-SM because dp will find the *lowest* total cost
 [p,q,C] = dp(1-SM);
 
 % Overlay the path on the local similarity matrix
 %hold on; plot(q,p,'r'); hold off
 
 % Path visibly follows the dark stripe
 
 % Plot the minimum-cost-to-this point matrix too
 %subplot(122)
 %imagesc(C)
 %hold on; plot(q,p,'r'); hold off
 
 % Bottom right corner of C gives cost of minimum-cost alignment of the two
 C(size(C,1),size(C,2));

 % This is the value we would compare between different 
 % templates if we were doing classification.
 
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
 d2x = resize(d2x', length(d1),1);
 
end