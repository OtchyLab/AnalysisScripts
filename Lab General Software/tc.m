for k = 1:(size(J34_syllable_tc,2))
    j34_timecourse(k)=J34_syllable_tc(1,k)*24*60+J34_syllable_tc(2,k)*60+J34_syllable_tc(3,k);
end;
nbins=abs(j34_timecourse(size(J34_syllable_tc,2))-j34_timecourse(1));
n=hist(j34_timecourse,floor(nbins/10));
figure; plot(n);