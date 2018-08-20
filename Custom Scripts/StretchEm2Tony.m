%Script converts filtered audio and neuro data matrices to wav files for
%Tony to process.

PPaudio;
PPneuro;

filLoc = 'C:\Users\Tim\Desktop\DirUndir XCorr Sets\For Fantana Figure\Global Xcorr Matrix\TonyData';
prefix = 'PUR692_04_20_';

if length(PPaudio)~=length(PPneuro)
    print('Cell arrays should be the same size. Double check the inputs');
    return
end

numRend = length(PPaudio);

for j=1:numRend
   nameAUD = [prefix '_r' num2str(j) '_AUD.wav'];
   wavAud = PPaudio{j}./max(abs(PPaudio{j}));
   wavwrite(wavAud,44150,[filLoc '\' nameAUD]);
   
%    nameNEU = [prefix '_r' num2str(j) '_NEU.wav'];
%    wavNeu = PPneuro{j}./max(abs(PPneuro{j}));
%    wavwrite(wavNeu,44150,[filLoc '\' nameNEU]);
end
