function [a,b] = golay(n);
%
% create a pair of golay codes of order n
% use assorted randomizations which preserve 
% complementary series properties
%
%IN: n = order of golay sequences
%
%OUT: a = golay sequence 1
%     b = golay sequence 2
%
%note : xcorr(a) = xcorr(b) = delta_fn
%
% AL, Caltech, 3/01
%

a1 = sign(unifrnd(-1,1,1,1));
a2 = sign(unifrnd(-1,1,1,1));

a = [a1 a2];
b = [a1 -a2];

for i = 2:n
	a2 = [a b];
	b = sign(unifrnd(-1,1,1,1))*[a -b];
	a = sign(unifrnd(-1,1,1,1))*a2;
	if sign(unifrnd(-1,1,1,1)) > 0
		a = fliplr(a);
	end;
	if sign(unifrnd(-1,1,1,1)) > 0
		b = fliplr(b);
	end;
end;

