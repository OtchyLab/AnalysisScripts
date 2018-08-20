function r = findMotifsInAudio(audio, templateAudio, sampRate)

%boxFiltLen = round(.010*sampRate);
%boxFilter = repmat(1/boxFiltLen, boxFiltLen,1);

audio = audio - mean(audio);
%audio = audio.^2;
%audio = conv(audio, boxFilter);

templateAudio = templateAudio - mean(templateAudio);
%templateAudio = templateAudio.^2;
%templateAudio = conv(templateAudio, boxFilter);

[r, lags] = xcorr(audio, templateAudio);
%r = r(ceil(length(r/2)):end);
%find(r>
%figure(2); plot(audio);
%figure(3); plot(templateAudio);
figure(4); plot(lags/40000, r);