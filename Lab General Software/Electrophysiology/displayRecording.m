function displayRecording(exper,num)

subplot(length(exper.sigCh)+1, 1, 1);
audio = loadAudio(exper,num);
specgram(audio);
for(ch = 1:length(exper.sigCh))
    subplot(length(exper.sigCh)+1, 1, ch+1);
    data = loadData(exper,num,exper.sigCh(ch));
    plot(data);
end
