function SCRAMBLED_SIGNAL=markolab_phase_scramble_1d(SIGNAL,DEBUG)
%phase scramble a 1D signal
%
%

% take the fft

sig_fft=fft(SIGNAL);

% get amplitude

sig_amp=abs(sig_fft);

% get theta

sig_theta=angle(sig_fft);

% scrambled theta

scr_theta=angle(fft(rand(size(SIGNAL))));

% new signal

SCRAMBLED_SIGNAL=real(ifft(sig_amp.*exp(i.*scr_theta)));

if DEBUG

	% original signal

	figure();
	ax(1)=subplot(4,1,1);
	original_signal=real(ifft(sig_amp.*exp(i.*sig_theta)));
	plot(original_signal);
	title('Original signal');

	% scrambled signal

	ax(2)=subplot(4,1,2);
	plot(SCRAMBLED_SIGNAL);
	title('Phase scrambled signal');

	% two-sided spectrum of original signal

	ax(3)=subplot(4,1,3);
	plot(sig_amp);
	title('abs(FFT) original signal');

	% two-sided spectrum of scrambled signal

	ax(4)=subplot(4,1,4);
	plot(abs(fft(SCRAMBLED_SIGNAL)));
	title('abs(FFT) phase scrambled signal');

	linkaxes(ax,'x');
	axis tight

end
