function [pitch pitchGoodness] = sparse_fftm(data,centerFreq,numHarmonics,fs)

%Inputs
% data: an n x 1 column vector of the waveform

% centerFreq: an estimate of roughly where the pitch would be; it does not
% have to be exact but +/- 100Hz (e.g., 520Hz)

% numHarmonics: the number of harmonics including the fundamental (pitch) to use to compute pitch;
% typically, 20 is used (e.g., 520Hz, 1040Hz, 1560Hz, etc.)

% fs: sampling frequency

% Outputs
% pitch: best estimate of pitch (Hz)

% pitchGoodness: how much of total power in the data the pitch its
% harmonics can explain; 1 would mean the data is a pure harmonic stack


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%calculate waveform duration and other necessary parameters for subsequent
%fft
len = length(data);
sampleFreqRes = 1./(len./fs);

centerFreqidx = round(centerFreq./sampleFreqRes);
changeResCenter = (sampleFreqRes./len).*(centerFreqidx);
numPointsTest = round(sampleFreqRes./changeResCenter);

centerFreqstep = centerFreq./centerFreqidx;

temp2 = (sampleFreqRes./changeResCenter)./2;
temp1=round(((1./(centerFreq./centerFreqidx)).*fs)+temp2);

if temp1 > len
centerFreqidx = centerFreqidx - 1;    
centerFreqstep = centerFreq./centerFreqidx;
numPoints=round((1./(centerFreqstep).*fs)+temp2);

else
 
numPoints =   ((1./(centerFreq./centerFreqidx)).*fs);  
    
end
clearvars temp*

data_in = data(1:numPoints);   

%main step of fft and pitch estimation using sparse assumption
for i=1:numPointsTest
    
 data_temp = data_in(i:end);   
 
 temp1{i,1} = abs(fft(data_temp));
 temp2(i,1) = sum(abs(fft(data_temp)));
 for j = 1:numHarmonics
 temp_idx  = centerFreqidx.*(j);  
 temp3{i,1}(j,1)= temp1{i,1}(temp_idx+1);
     
 end
 
 temp4(i,1)=sum(temp3{i,1});
 temp5(i,1)=temp4(i,1)./temp2(i,1);
 
 
end
[temp6,temp_idx] = max(temp5);

pitchGoodness = temp6.*2;
pitch = ((fs./(numPoints - (temp_idx-1)))).*centerFreqidx;

end