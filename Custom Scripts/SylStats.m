function [m,s] = SylStats(data)

mask = (data(:,4)>0 & data(:,4)<100);
m = mean(data(mask,4));
s = std(data(mask,4));

% [m(1),s(1)]=SylStats(data)
% [m(length(m)+1),s(length(s)+1)]=SylStats(data)
% DayM = mean(m)
% DayStdM = std(m)

end