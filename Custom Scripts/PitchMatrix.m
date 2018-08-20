function [pitchMat] = PitchMatrix
%This script loads the dataset from the saved process file (from StretchEm)
%and returns:
%pitchMat = m x n matrix of all syllable pitches
%       m = renditions
%       n = syllable number
%keys = the filesnames of each rendition
%startSyl = the index of the starting syllable of sequence for the file
%sequence = the syllable sequence represented in the matrix

%Request and load dataset
[fname,pathLoc] = uigetfile('C:\Users\Tim\Desktop\pCAF Deafen data\Pur822\*.mat');
tots = [pathLoc, fname];
load(tots,'data','filenames','sequence')

%Select target syllable bounds
a = transpose(data.templatesyllBreaks(:,1));
%a = a(:);

%Syllable targeting
syls = [2,4,5]; %Which syllables to measure pitch
pOffset = round((20/1000)*44150); %Distance from syllable onset (in samples) to measure pitch
pLength = 220; %Length (in samples) of syllable segment to analyze

%Pitch Calc Parameters
centerFreq = 570;
numHarmonics = 20;
fs = 44150;

%Interval Durations from the paths and template files
rendNum = length(data.audio);
pitchMat = [];
for i = 1:rendNum
    %Rover rendition path and rendition length
    path = [data.p{i},data.q{i}];
    rendLength = size(data.audio{i},2);
    
    %Calculate pitch and add to the stack
    [warpedOut] = getWarpedStarts(path,a(:));
    
    %Step through each syllable and calc pitch, or give NaN
    for j = 1:size(data.templatesyllBreaks,1)
        if ~isempty(syls(j==syls)) %if it's one of the chosen
            %Get audio snippet
            startPnt = round(fs*(warpedOut(j)/1000))+pOffset;
            indx = startPnt:(startPnt+pLength-1);
            snip = data.audio{i}(indx);
            %Calculate pitch
            [pitchRend(j) pitchGoodness] = sparse_fftm(snip,centerFreq,numHarmonics,fs);
        else
            %Fill with dummy variable
            pitchRend(j) = NaN;
        end
    end
    pitchMat = [pitchMat; pitchRend];
end

%Plot interval output stats
pitch.m = mean(pitchMat,1);
pitch.std = std(pitchMat);
pitch.cv = pitch.std./pitch.m;

fig1 = figure;
subplot(2,1,1)
hold on
scatter(1:length(pitch.m),pitch.m,'or');
errorbar(1:length(pitch.m),pitch.m,pitch.std,'or');
set(gca,'XTick',[],'XTickLabel',[])
ylabel('Syllable Pitch (Hz)')
xlim([0,length(pitch.m)+1])
tits = ['Syllable Pitch for ' fname(1:end-11)];
title(tits,'Interpreter','none')

subplot(2,1,2)
hold on
scatter(1:length(pitch.cv),pitch.cv,'sk')
line([0,length(pitch.cv)+1],[nanmean(pitch.cv(1:length(pitch.cv))),nanmean(pitch.cv(1:length(pitch.cv)))],'Color','k','LineStyle','--')
%line([0,length(pitch.cv)+1],[mean(pitch.cv(2:2:length(pitch.cv))),mean(pitch.cv(2:2:length(pitch.cv)))],'Color','k','LineStyle','--')
temp = [];
nums = 1:.5:5;
for i = 1:2:length(pitch.cv)
    temp = [temp; char(num2str(nums(i)))];
    %temp = [temp; char('G')];
end
set(gca,'XTick',1:length(pitch.cv))%,'XTickLabel',temp)
xlabel('Syllable')
ylabel('CV')
xlim([0,length(pitch.m)+1])

%Save data to file
savename = [pathLoc, 'pitch\' fname(1:end-11), ' pitch.mat'];
save(savename,'pitchMat','filenames','sequence','pitch')

saveas(fig1,[pathLoc, 'pitch\' fname(1:end-11), '.fig'])
%saveas(fig1,[pathLoc, 'pitch\' fname(1:end-11), '.eps'])