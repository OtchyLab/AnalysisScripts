function [LinPath,base_syllBreaks,future_syllBreaks] = HVCwarp(baseline)%,future,template,templatesyllBreaks)
%baseline - mean HVC power trace from the baseline case
%future - mean HVC power trace for some future day that will be compared to
%   the baseline
%template - the spectrogram/template for the baseline
%templatesyllBreaks - the timebins that demarcate the syllable boundaries
%   in the baseline

%Get baseline power trace
% [fname loc] = uigetfile('C:\Users\Tim\Desktop\For Ali Paper Test Alignments\Pur639\','Select baseline dataset source');
% load([loc fname],'mNeuro')
% baseline = mNeuro;

%Get future power trace
[fname loc] = uigetfile('C:\Users\Tim\Desktop\For Ali Paper Test Alignments\Pur639\','Select future dataset source');
load([loc fname],'mNeuro')
future = (smooth(mNeuro,220))';

%Get template info
[fname loc] = uigetfile('C:\Users\Tim\Desktop\For Ali Paper Test Alignments\Pur639\','Select baseline template source source');
load([loc fname])
template = template;
templatesyllBreaks = templatesyllBreaks;

% This will warp the HVC traces from the future day to the baseline day
[val, p] = DTWsubBand(baseline,future);
futureRealign = alignSeries(future,p);

%Create anchor point vectors -- in this case syl/gap edges
base_syllBreaks = [];
for ind = 1:size(templatesyllBreaks,1)     
     base_syllBreaks = [base_syllBreaks, templatesyllBreaks(ind,:)];
end
base_syllBreaks = [1, base_syllBreaks, size(template,2)];

future_syllBreaks = getRemapAnchors(p,44*base_syllBreaks)./44; %correct for bin size mismatch

%Create the full linear path
LinPath = [];
for j = 2:length(base_syllBreaks)
    path_t = linspace(future_syllBreaks(j-1),future_syllBreaks(j),base_syllBreaks(j)-base_syllBreaks(j-1)+1);
    if j~=length(base_syllBreaks)
        LinPath = [LinPath,path_t(1:end-1)];
    else
        LinPath = [LinPath,path_t];
    end
end

%Plot output
figure(66)
subplot(2,2,1)
hold on
plot(LinPath-(1:length(LinPath)))
line([0,length(LinPath)],[0,0],'Color','y');

subplot(2,2,2)
hold on
plot(diff(LinPath-(1:length(LinPath))))

subplot(2,2,3:4)
hold on
plot(baseline)
plot(future,'r')
plot(futureRealign,'g')




