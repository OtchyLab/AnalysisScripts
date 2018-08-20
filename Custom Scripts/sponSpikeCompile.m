function sponSpikeCompile

%List of spike files to load
% spikeFileList = [{'Grn046_140515_spon2_spikeRate.mat'},...   %Grn046
%                         {'Grn046_140515_spon3_spikeRate.mat'},...
%                         {'Grn046_140515_spon4_spikeRate.mat'}];
% lesionTOD = 11.5; 
% lightsOut = 22;
% lightsOn = 33;
% lightsOut2 = 46;

spikeFileList = [{'Grn121_140906_spon2_spikeRate.mat'},...   %Grn121
                        {'Grn121_140906_spon3_spikeRate.mat'},...
                        {'Grn121_140906_spon4_spikeRate.mat'}];
lesionTOD = 14.25; 
lightsOut = 23.;
lightsOn = 34;
lightsOut2 = 46;

% spikeFileList = [{'Grn141_141210_spon2_spikeRate.mat'},...   %Grn141
%                         {'Grn141_141210_spon3_spikeRate.mat'},...
%                         {'Grn141_141210_spon4_spikeRate.mat'}];
% lesionTOD = 12; 
% lightsOut = 22;
% lightsOn = 35;
% lightsOut2 = 44;

% spikeFileList = [{'Grn186_1500424_spon2_spikeRate.mat'},...   %Grn186
%                         {'Grn186_1500424_spon3_spikeRate.mat'},...
%                         {'Grn186_1500424_spon4_spikeRate.mat'}];
% lesionTOD = 11.5; 
% lightsOut = 23.5;
% lightsOn = 34;

sumRates = [];
interpRates = [];
for i = 1:length(spikeFileList)
    %Load in data
    load(spikeFileList{i}, 'fileHour', 'spikeRate')
    fileHour = fileHour;
    sumRates(i,:) = spikeRate;
    
    %Zero time at lesion time
    relTime = fileHour - lesionTOD;
    
    %Interpolate between points for later recombining
    interpTimes = floor(min(relTime)):0.25:ceil(max(relTime));
    interpRates(i,:) = interp1(relTime, spikeRate, interpTimes,'linear','extrap');
end

%Fix numbering wrt to lesion time
relLightsOut = lightsOut - lesionTOD;
relLightsOn = lightsOn - lesionTOD;
relLightsOut2 = lightsOut2 - lesionTOD;

%Generate parsing indices
preIdx = interpTimes < 0;
postIdx = (interpTimes > 0) & (interpTimes < relLightsOut);
post1Idx = (interpTimes > relLightsOn) & (interpTimes < relLightsOut2);

%Normalize by the mean of the prelesion values
normRates = [];
for i = 1:length(spikeFileList)
    norm = mean(interpRates(i,preIdx));
    normRates(i,:) = interpRates(i,:)./norm;
end

%Take average across all channels
meanRate = mean(normRates,1);

%Plot the output
figure(10); clf
plot(interpTimes, normRates(1,:), '-xk'); hold on
plot(interpTimes, normRates(2,:), '-xb');
plot(interpTimes, normRates(3,:), '-xr');
plot(interpTimes, meanRate, '-sg', 'LineWidth', 2);
axis tight; y = ylim;
line([0, 0], [0, y(2)], 'Color', 'r') %Lesion time
line([relLightsOut, relLightsOut], [0, y(2)], 'Color', 'k') %Lights out
line([relLightsOn, relLightsOn], [0, y(2)], 'Color', 'y') % Lights on
line([relLightsOut2, relLightsOut2], [0, y(2)], 'Color', 'k') %Lights out


%Parse data
preIdx = interpTimes < 0;
postIdx = (interpTimes > 0) & (interpTimes < relLightsOut);
post1Idx = (interpTimes > relLightsOn) & (interpTimes < relLightsOut2);

preTime = interpTimes(preIdx);
preRate = meanRate(preIdx);

postTime = interpTimes(postIdx);
postRate = meanRate(postIdx);

post1Time = interpTimes(post1Idx);
post1Rate = meanRate(post1Idx);


%Saving to file
output.Name = spikeFileList;
output.preTime = preTime;
output.preRate = preRate;
output.postTime = postTime;
output.postRate = postRate;
output.post1Time = post1Time;
output.post1Rate = post1Rate;

%%%%%%%%%%%%%%%%%%%%%
% Save to output location
%%%%%%%%%%%%%%%%%%%%%
%Function output location
saveLoc = 'C:\Users\Tim\Desktop\Nif Project Figures\ElectroLesions\Spontaneous Analysis\Summary Spontaneous Data 05202015.mat';
m = exist(saveLoc);
if m ==2 
    %File already exists
    clear('sponSumStats');
    load(saveLoc, 'sponSumStats')
    sponSumStats(end+1) = output;
else
    %No file yet created
    sponSumStats = output;
end

%Save the updated data to file
save(saveLoc, 'sponSumStats')

%Write the noisy file names to a text file for removal later
display('done')
