%% Select Data
clear all
%Process one or more DAT files by name
[temp,path] = uigetfile('*.mat','Select the interval files to process','MultiSelect','on');
cd(path);
if size(temp,2)>1
    for i = 1:size(temp,2)
        fnames(i).name = char(temp{i});
    end
else
    fnames.name = temp;
end
clear('temp')
numFiles = length(fnames);

%% Load and process data sequentially
means = [];
stds = [];
n_means = [];
n_stds = [];
cvs = [];
for i = 1:numFiles
    load(fnames(i).name,'IntMat','intervals')
    if i == 1
        numInts = size(IntMat,2);
    else
        if numInts ~= size(IntMat,2)
            print(['Number of intervals in file ' fnames(i).name ' does not match first file interval number']);
            return
        end
    end
    all_IntMat{i} = IntMat;
    all_Intervals(i) = intervals;
    clear('IntMat','intervals')
    
    %Reorg data for plotting ease
    means = [means,all_Intervals(i).m'];
    n_means = [n_means,(all_Intervals(i).m./all_Intervals(1).m)'];
    stds = [stds,all_Intervals(i).std'];
    n_stds = [n_stds,(all_Intervals(i).std./all_Intervals(1).m)'];
    cvs = [cvs,all_Intervals(i).cv'];
    
    %Create elapsed time array
    splits = regexp(fnames(i).name,'_','split');
    numdate = datenum(splits(2),'yymmdd');
    if i == 1
        day(i) = 1;
        basedate = numdate;
    else
        day(i) = numdate-basedate+1;
    end
end

%% Plot data
%Set the intervals to plot
ints = 1:numInts;
%ints = 5:6;
labels = {'1','1G2','2','2G3','3','3G4','4','4G5','5','5G6','6'};
cols = {'b','r','g','k','m','c','y'};
jams = [];

figure;
%Plot real interval durations
subplot(3,1,1)
hold on
%plot(day,means(ints,:)')
for i = ints
    errorbar(day,means(i,:)',stds(i,:)',cols{i})
end
xlim([min(day),max(day)])
ylims = ylim;
for i = 1:length(jams)
    line([jams(i),jams(i)],ylims,'LineStyle',':')
end
ylabel('Interval Duration (ms)')
set(gca,'XTick',day,'TickDir','out');
set(gca,'Box','off')
tits = ['Intervals for ' fnames(1).name(1:end-10) ' through ' fnames(end).name(1:end-10)];
title(tits,'Interpreter','none')

%Plot normalized interval durations
subplot(3,1,2)
hold on
%plot(day,n_means(ints,:)')
for i = ints
    errorbar(day,n_means(i,:)',n_stds(i,:)',cols{i})
end
xlim([min(day),max(day)])
ylims = ylim;
for i = 1:length(jams)
    line([jams(i),jams(i)],ylims,'LineStyle',':')
end
xlim([min(day),max(day)])
ylabel('Interval Duration (norm)')
set(gca,'XTick',day,'TickDir','out');
set(gca,'Box','off')

subplot(3,1,3)
hold on
for i = ints
    plot(day,cvs(i,:)',cols{i},'DisplayName',labels(i))
end
xlim([min(day),max(day)])
ylims = ylim;
for i = 1:length(jams)
    line([jams(i),jams(i)],ylims,'LineStyle',':')
end
ylabel('CV')
xlabel('Days')
set(gca,'XTick',day,'TickDir','out');
set(gca,'Box','off')

