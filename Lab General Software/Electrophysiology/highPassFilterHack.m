function cleanData = highPassFilterHack(data, targetFreqHz, sampleRate)

%freq is measured as number of periods within the length of data.
fftdata = fft(data(:,1),length(data));
targetFreq = floor((length(data)/sampleRate)*targetFreqHz); %targetFreq = #ofSeconds in data times 60 cycles per second.
targetRolloff = floor((length(data)/sampleRate)*500)
filter = zeros(length(fftdata), 1);
filter(targetFreq+1:targetFreq+targetRolloff+1) = sin(([0:targetRolloff]/targetRolloff)*pi/2).^2;
filter(length(data)-(targetFreq+targetRolloff-1):length(data)-(targetFreq-1)) = sin(([targetRolloff:-1:0]/targetRolloff)*pi/2).^2;
filter(targetFreq+targetRolloff+1+1:length(data)-(targetFreq+targetRolloff-1)-1) = 1;
fftdata = fftdata .* filter;

%fftdata(rangeLow+1:rangeHigh+1) = 0; %indice of freq is 1+freq.
%fftdata(length(data)-(rangeHigh-1):length(data)-(rangeLow-1)) = 0;
cleanData = ifft(fftdata);
