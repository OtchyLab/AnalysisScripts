function fraction = GetFrac(spikeI)

burst =0;

for j=1:length(spikeI)-1
    if (spikeI(j)>0.0067 && spikeI(j+1)>(0.0067))
            burst=burst+1;
    end
end
    fraction=1-(burst/length(spikeI));
end