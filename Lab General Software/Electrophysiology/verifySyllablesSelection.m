function [bValid, syllableStartNdx, syllableEndNdx] = verifySyllablesSelection(audio, sampRate, syllableStartNdx, syllableEndNdx)

h = figure(1000);
displayAudioSpecgram(audio, sampRate, 0);
hold on;
startTime = (syllableStartNdx - 1)/sampRate;
endTime = (syllableEndNdx - 1)/ sampRate;
for(nSyll = 1:length(startTime))
    lines(nSyll, 1) = line([startTime(nSyll), startTime(nSyll)], ylim, 'Color', 'yellow');
    lines(nSyll, 2) = line([endTime(nSyll), endTime(nSyll)], ylim, 'Color', 'red');
end

char = input('Are syllables markers correct?', 's');
bValid = (char == 'y' | char == 'Y');

close(h);

