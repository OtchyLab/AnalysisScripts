%Process descriptive stats (maybe figures) for the nanoclips echem data

%Reset workspace
clear

%Load the raw data
load('/Users/Tim/Desktop/PadEchem.mat')

%Initialize output structure
stats = [];

[numPads, numVars] = size(PadEchem);

%Process descriptive stats sequentially
for i = 1:numVars
    %For convenience, convert from tabel structure
    t = table2array(PadEchem(:,i));
    
    %Strip out the NaNs
    t = t(~isnan(t));
    numel(t)
    
    %Mean
    stats.mean(i) = mean(t);
    
    %Median
    stats.median(i) = median(t);
    
    %Std
    stats.std(i) = std(t);
    
    %Mean
    stats.sem(i) = std(t)/sqrt(numel(t));
    
    %Range
    stats.range(i,:) = [min(t), max(t)];
   
    %Copy out
    v{i} = t;
end

%Plot output
figure(69); clf

%Scatter plot the Z (log scale)
subplot(1,2,1)

%Plot preZ
jit = (rand(size(v{1}))-0.5)./5;
z1 = scatter(1.*ones(size(v{1}))+jit, v{1}.*1000); hold on

jit = (rand(size(v{3}))-0.5)./5;
z2 = scatter(2.*ones(size(v{3}))+jit, v{3}.*1000);

xlim([0.5,2.5])
ylabel('Impedance (\Omega)')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1,2], 'XTickLabel', {'Bare Gold'; 'EIROF'})
set(gca,'yscale','log')


%Scatter plot the CSC
subplot(1,2,2)

%Plot preZ
jit = (rand(size(v{2}))-0.5)./5;
csc1 = scatter(1.*ones(size(v{2}))+jit, v{2}.*.8); hold on

jit = (rand(size(v{4}))-0.5)./5;
csc2 = scatter(2.*ones(size(v{4}))+jit, v{4});

xlim([0.5,2.5])
% ylim([10e-1, 10e1])
ylabel('CSC (\muC/cm^2)')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1,2], 'XTickLabel', {'Bare Gold'; 'EIROF'}, 'YTick', [0, 50])

a = 1;



