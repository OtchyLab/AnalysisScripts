function a=autocorr1(spiketrain, resolution, length);

diff=[];
time=max(spiketrain)-min(spiketrain);
fr=size(spiketrain,1)/time;
for i=1:size(spiketrain,1)
    for j=i+1:size(spiketrain,1)
    diff(end+1)=spiketrain(j)-spiketrain(i);
    end
end
x=(0:resolution:length);
if (isempty(diff))
    diff(1:size(x,2))=0;
end

