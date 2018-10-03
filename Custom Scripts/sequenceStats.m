%function sequenceStats
%Function will calculate and plot a few different representations of
%sequence variability. 
%
%The main input folder is selected by the user using a GUI prompt. in that
%folder should be several saved sequence matrices created by StretchEm or
%StretchEmLite. Each file in that folder should be named:
%     BirdName_YYMMDD_SequenceStats.mat
%
% Created 9/28/2018 TMO

%% Load from file

%Select the data folder to process
parent = uigetdir(cd, 'Which folder to process?');

%Get the list of files to process
seqFiles = dir([parent, filesep, '*SequenceStats.mat']);

%Figure out how many files you're going to process
numFiles = numel(seqFiles);

%Loop through the list of files and copy data from the harddrive to the
%local variable
seqs = [];
for i = 1:numFiles
    %Select file to load
    filename = seqFiles(i).name;

    %Load from the file
    load([parent, filesep, filename], 'seqMatrix');
    
    %Copy to local variable
    seqs(i,:,:) = seqMatrix;
    
    %Delete the loaded variable
    clear('seqMatrix')
end

%% Extract date information for later plotting
surgeryDate = '180111'; %change this as necessary
surgeryNum = datenum(surgeryDate, 'yymmdd');
days = [];
for i = 1:numFiles
    %Select file to load
    filename = seqFiles(i).name;

    %Split up the filename using the underscore
    sp = regexp(filename, '_', 'split');
    
    %Select the portion that includes the date
    dStr = sp{2};
    
    %Convert to a datenum for ease
    dNum = datenum(dStr, 'yymmdd');
    
    %Day before or since surgery
    days(i) = dNum - surgeryNum;

end

%% Simple shannon entropy

%Loop through the each day of data and calculate simple shannon entropy
shanH = [];
for i = 1:numFiles
    %Select out the day of data to work with
    t = squeeze(seqs(i,:,:));
    
    %Sum the number of transitions in the whole matrix
    transNum = sum(t(:));

    %Convert count matrix into probability distribution
    t_pd = t./transNum;

    %Fix the log(0) problem
    t_pd(t_pd == 0) = eps;
    
    %Calculate entropy for each sequence outcome
    t_H = -t_pd.*log2(t_pd);
    
    %Sum over transitions
    shanH(i) = sum(t_H(:));
    
end

%Calculate maximum entropy
[m, n] = size(t_H);
maxH = log2(m*n);

%Plot the result

%Initialize figure
figure(100); clf

%Plot the timeseries
q1 = plot(days, shanH, '-o');
hold on %prevent overwrite

%Format the ts plot
q1.LineWidth = 1.5;

%Plot the max entropy
q2 = line([days(1), days(end)], [maxH, maxH]);

%Format the max line
q2.Color = 'k';
q2.LineStyle = '--';

%Plot surgery line
q3 = line([0,0], [0,5]);

%Format the surgery line
q3.Color = 'r';
q3.LineStyle = ':';

%Format the figure
xlim([-5, 35])
set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', 0:5, 'LineWidth', 1.5) %fix ugly axes
xlabel('Time (days)')
ylabel('Entropy (bits)')
title('Shannon Entropy of Sequences')

%Set the size and location of the figure
set(gcf, 'Units', 'Inches', 'Position', [4, 4, 4, 4])


%% Node Transition Probability Plots

%Loop through the each day of data and calculate transition probabilities
%from each node (i.e., syllable)
[m, n] = size(squeeze(seqs(1,:,:)));
sylPrev = [];
Pij = [];
for i = 1:numFiles
    %Select out the day of data to work with
    t = squeeze(seqs(i,:,:));
    
    %Sum the number of transitions in the whole matrix
    totalNum = sum(t(:));
    
    %Sum the number of transitions for each node
    transNum = sum(t, 2);

    %Convert count matrix into trans probabilities for each node
    t_pd = t./transNum;
    
    %Syllable prevalence
    sylPrev(i,:) = transNum./totalNum;

    %Dump to running variable
    Pij(i,:,:) = t_pd;
    
end

%Calculate normalized to the mean of the first 3 days (i.e., pre-surgery)
d13 = sylPrev(1:3,:);
m_d13 = mean(d13,1);
normPrev = sylPrev./m_d13;

%Plot the results

%Initialize syllable prevalence figure
figure(101); clf

%Plot the prevalence for each syllable
subplot(1,2,1)
r = plot(days, sylPrev, '-o', 'Linewidth', 1.5);

%Plot surgery line
q = line([0,0], [0,1]);

%Format the surgery line
q.Color = 'r';
q.LineStyle = ':';

%Format the figure
xlim([-5, 35])
set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', 0:0.5:1, 'LineWidth', 1.5) %fix ugly axes
xlabel('Time (days)')
ylabel(['Syllable Prevalence (%)'])

%Add a legend
lgd = legend('Syl 1','2','3', '4', '5');

%Plot the normalized prevalence
subplot(1,2,2)
r = plot(days, normPrev, '-o', 'Linewidth', 1.5);

%Plot surgery line
q = line([0,0], [0,2]);

%Format the surgery line
q.Color = 'r';
q.LineStyle = ':';

%Format the figure
xlim([-5, 35])
set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', 0:2, 'LineWidth', 1.5) %fix ugly axes
xlabel('Time (days)')
ylabel(['Norm Prevalence'])

%Set the size and location of the figure
set(gcf, 'Units', 'Inches', 'Position', [4.5, 10, 7.25, 3])

%Initialize node trans figure
figure(102); clf

%Plot the timeseries for each node transition probability
for i = 1:m
    %Select the figure subplot
    subplot(m, 1, i)
    
    %Select out the node to plot
    u = squeeze(Pij(:,i,:));
    
    %Plot the probabilities for the node
    r = plot(days, u, '-o', 'Linewidth', 1.5);
    
    %Plot surgery line
    q = line([0,0], [0,1]);
    
    %Format the surgery line
    q.Color = 'r';
    q.LineStyle = ':';
    
    %Format the figure
    xlim([-5, 35])
    set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', 0:0.5:1, 'LineWidth', 1.5) %fix ugly axes
    ylabel(['Trans from Syl ' num2str(i)])
    
end
xlabel('Time (days)')
lgd = legend('Syl 1','2','3', '4', '5');

%Set the size and location of the figure
set(gcf, 'Units', 'Inches', 'Position', [1, 6.5, 3.2, 7.5])

%% Test the independence of sequence nodes
% 
% %Loop through the each day of data and calculate 1st and 2nd order independence
% indep = [];
% for i = 1:numFiles
%     %Select out the day of data to work with
%     t = squeeze(seqs(i,:,:));
%     
%     %Sum the number of transitions in the whole matrix
%     transNum = sum(t(:));
% 
%     %Convert count matrix into probability distribution
%     t_pd = t./transNum;
% 
%     %Fix the log(0) problem
%     t_pd(t_pd == 0) = eps;
%     
%     %Calculate entropy for each sequence outcome
%     t_H = -t_pd.*log2(t_pd);
%     
%     %Sum over transitions
%     shanH(i) = sum(t_H(:));
%     
% end
% 
% %Calculate maximum entropy
% [m, n] = size(t_H);
% maxH = log2(m*n);
% 
% %Plot the result
% 
% %Initialize figure
% figure(100); clf
% 
% %Plot the timeseries
% q1 = plot(days, shanH, '-o');
% hold on %prevent overwrite
% 
% %Format the ts plot
% q1.LineWidth = 1.5;
% 
% %Plot the max entropy
% q2 = line([days(1), days(end)], [maxH, maxH]);
% 
% %Format the max line
% q2.Color = 'k';
% q2.LineStyle = '--';
% 
% %Plot surgery line
% q3 = line([0,0], [0,5]);
% 
% %Format the surgery line
% q3.Color = 'r';
% q3.LineStyle = ':';
% 
% %Format the figure
% xlim([-5, 35])
% set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', 0:5, 'LineWidth', 1.5) %fix ugly axes
% xlabel('Time (days)')
% ylabel('Entropy (bits)')
% title('Shannon Entropy of Sequences')
% 
% %Set the size and location of the figure
% set(gcf, 'Units', 'Inches', 'Position', [4, 4, 4, 4])
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
