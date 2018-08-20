function [I,a1,b1,a2,b2] = get_golay(x, ich, och)
%
% get golay input/output sequences
% from observer data matrix
%
% note: we assume the lengths of a and b 
%       are equal, ie. x is of even length
%
%IN: x = observer data matrix
%    ich = channel number for input time series
%    och = channel number for output time series
%
%OUT: I = impulse response 
%     a1 = golay input sequence a
%     b1 = golay input sequence b
%     a2 = golay output sequence a
%     b2 = golay output sequence b
%
% AL, Caltech, 3/01
%

N = size(x,2);
M = N/2;
z=find(abs(x(ich,1:M)./(eps+max(x(ich,1:M))))>.75);
probelen = z(end);

a1 = x(ich,1:probelen);
a2 = x(och,1:M);
b1 = x(ich,M+1:M+probelen);
b2 = x(och,M+1:2*M);

I = gimpulse(a1,b1,a2,b2);
