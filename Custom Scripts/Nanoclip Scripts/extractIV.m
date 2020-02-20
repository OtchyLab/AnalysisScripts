%Script for extracting acute stim parameters
clear all

%Working folder (hooks)
%mother = 'V:\SongbirdData\LLR32';  timeS = 153500; timeE = 170000;
%mother = 'V:\SongbirdData\LLR20';  timeS = 143000; timeE = 153000;

%Working folder (nc)
%mother = 'V:\SongbirdData\LLR32'; timeS = 170213; timeE = 173000;
mother = 'V:\SongbirdData\Rd40'; timeS = 153300; timeE = 172945;

%Load list of files from mother dir
files = dir([mother, filesep, 'stim*.mat']);

%Parse out the files that you actually care to analyze; toss the rest
filenames = []; times = [];
for i = 1:numel(files)
    x = files(i).name; %get names
    a = regexp(x, '_', 'split');
    filenames{i} = x;
    times(i) = str2num(a{3}(1:end-4));
end

mask = times > timeS & times < timeE;
strpFiles = filenames(mask);

%Sequenctially extract the data we care about
cur = []; vHigh = []; vLow = []; cHigh = []; cLow = [];
for i = 1:numel(strpFiles)
    load([mother, filesep, strpFiles{i}]);

    %Extract mean current and voltagre
    cur(i) = data.stim.current_uA;
%     vLow(i) = data.voltage_range(1); %original
%     vHigh(i) = data.voltage_range(2); %original

    %extract trial by trial current and voltagre
    for j = 1:size(data.ni.stim,1)
        %Voltage
        vHigh = [vHigh, max(data.ni.stim(j,:,1))];
        vLow = [vLow, min(data.ni.stim(j,:,1))];

        %Current
        x = data.ni.stim(j,:,2);
        highC = x(x>1);
        lowC = x(x<-1);
        cHigh = [cHigh, median(highC)];
        cLow = [cLow, median(lowC)];
        
    end
    
    clear('data')
end

%Organize that data in some useful way
curAll = [cLow, cHigh];
vAll = [vLow, vHigh];

%Threshold outliers (noise/errors)
vMask1 = vAll < 1 & vAll > -0.6; %toss
vMask2 = vAll > 6 | vAll < -6; %toss
keepMask = ~(vMask1 | vMask2);
curPlotA = curAll(keepMask);
vPlotA = vAll(keepMask);

curPlot = curPlotA(~selTot');
vPlot = vPlotA(~selTot');

%Plot that data
figure(27); clf
plot(curPlot, vPlot, 'ok')
xlim([floor(1.1*unique(min(curPlot))), ceil(1.1*unique(max(curPlot)))])
set(gca, 'Box', 'off', 'TickDir', 'out')
xlabel('Current (uA)')
ylabel('Voltage (V)')







