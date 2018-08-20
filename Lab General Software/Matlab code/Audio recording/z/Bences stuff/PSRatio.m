function y=PSRatio(scanrate,spectra,t,f);
%looks at the power in the spectrogram at two different bands and compares
%these


LowL=200; %lower band goes from LowL-LowH Hz
LowH=500;
HighL=2000;%High band goes from HighL-HighH Hz
HighH=5000;

deltaF=f(2)-f(1);
LL=ceil(LowL/deltaF); %index for LowL
LH=ceil(LowH/deltaF);
HL=ceil(HighL/deltaF);
HH=ceil(HighH/deltaF);
for k = 1:length(t)
    LowBand(k)=sum(spectra(LL:LH,k))/(LH-LL); 
    HighBand(k)=sum(spectra(HL:HH,k))/(HH-HL); 
end

y=HighBand./LowBand;
