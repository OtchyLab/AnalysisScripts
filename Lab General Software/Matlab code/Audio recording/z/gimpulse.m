function [I,H] = gimpulse(a,b,aout,bout)
%
% get the golay-probe based impulse response
%
%IN: a = input a
%    b = input b
%    aout = output a
%    bout = output b
%
%OUT: I = impulse response
%     H = transfer function
%
% AL, Caltech, 3/01
%
% [I,H] = gimpulse(a,b,aout,bout)
%

 
N = length(a);
M = length(aout);
zextra = 0; % powers of 2 extra zero-padding
nfft = 2^nextpow2(M+zextra);
ai = zeros(1,nfft);
bi = zeros(1,nfft);
ao = zeros(1,nfft);
bo = zeros(1,nfft);

% de-mean data
ai(1:N) = a - mean(a);
bi(1:N) = b - mean(b);
ao(1:M) = aout - mean(aout);
bo(1:M) = bout - mean(bout);

fai = fft(ai,nfft);
fbi = fft(bi,nfft);
fao = fft(ao,nfft);
fbo = fft(bo,nfft);

a = abs(a);
a = a(find(a > 0.1));
alen = length(a);
L2 = 2*alen*(mean(a)^2);
%L2 = mean(fai.*conj(fai) + fbi.*conj(fbi));

% divide by 2 also?
H = (fao.*conj(fai) + fbo.*conj(fbi))./L2;
I = real(ifft(H,nfft));
%H = log10(abs(H(1:round(length(H)/2)-1)));


