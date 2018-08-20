function [db,cancel] = cancelfeedback(I, tin, tout,show)

tin = tin - mean(tin);
tout = tout - mean(tout);
nlen = length(tout);
N = length(tout);
p1 = conv(tin,-I);
predict = zeros(1,nlen);
M = length(p1);
if M > nlen
   predict = p1(1:nlen);
else predict(1:M)  = p1;
end;
cancel = tout + predict;
db = 10*log10(sum(tout.^2)/sum(cancel.^2));

if show
	clf;
	plot(tout,'b-');
	hold on;
	plot(-predict,'r-');
	plot(cancel,'k');
	zoom on;
	sprintf('Signal power reduction of %d dB',db)
end;

