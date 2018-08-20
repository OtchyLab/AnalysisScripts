%function song=tcrossings(ratio,t,threshold);
%looks at the power in the spectrogram at two different bands and compares
%these
rate_threshold=3.5; %how large should the threshold be to call it a song?
s=find(ratio>threshold); 
d=diff(s);
deltaT=t(2)-t(1);


t_cross=find(d>1 & d<ceil(0.6/deltaT)); %find threshold crossings that are less than 0.6 seconds apart
if (length(t_cross)>2 & (s(t_cross(end)+1)-s(t_cross(1)))*deltaT>0.5) %song needs to consist of more than three syllables and be longer than 0.5 s
    silence=sum(d(1,[find(d>ceil(0.6/deltaT))])); %subtract long preiods of silence when calculating rate
    rate=length(t_cross)/((((s(end)-s(1))-silence))*deltaT);
    if (rate>rate_threshold & size(d(1,[find(d>ceil(0.6/deltaT))]))<2)
        song=1;
    else
        song=0;
    end
else
    song=0;
end
