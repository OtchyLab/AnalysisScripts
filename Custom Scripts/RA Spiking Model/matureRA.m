function [spikeStats, params] = matureRA(numHVC, inh_HVC_ratio, matureVecs)
% [spikeStats, params] = matureRA(numHVC, inh_HVC_ratio, matureVecs)
% written by TMO 06/11/2016; last modified 06/28/2016
%
% This script simulates the maturation of a model RA neuron (i.e., it
% models changes in spiking as a function of HVC and local inhibition.
% The code is based on the "RA_toy_model_lognormal_orig" script used in Garst-Orozco
% et al 2015.
%
% Inputs:
% numHVC:       scalar value for weight factor for HVC input  (def: 4)
% inh_HVC_ratio:scalar value for HVC-mediated inhibition (def: 200)
% matureVecs:   structure containing the vectures defining circuit maturation
%               includes:
%                   matureVecs.prune_vec
%                   matureVecs.mean_vec
%                   matureVecs.std_vec
%                   matureVecs.ages
%                   matureVecs.inh_vec
%                   matureVecs.HVC_vec
% Outputs:
% spikeStats:   structure of length that contains the trajectory of spiketrain statistics
%               includes:
%                   spikeStats.ages
%                   spikeStats.FR
%                   spikeStats.BF
%                   spikeStats.sparse
%                   spikeStats.spikeCorr
% params:       dunno
%

%Boolean switch for plotting the rasters
raster = true;

%These are (I think) the variables to search over
% inh_HVC_ratio = 200;
% numHVC = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%
%Constants for simulation
%%%%%%%%%%%%%%%%%%%%%%%%%
staticHVC = false; %Boolean switch for 
numTrials = 100; %Number of trials to run per age per neuron
N_HVC = 100; %Number of HVC neurons in complete pool
numLMAN = 2; %LMAN input multiplier
lineSize = 0.8; %raster plot line height
PSTH_bin_size = 20; %bin size in (ms)

%Simulation time vector
t_stop = 1000;  % length of each trial (ms);
dt = 0.2; %time step (ms)
t = dt:dt:t_stop;

%Setup the vectors defining the maturation of connectivity
prune_vec = matureVecs.prune_vec;
mean_vec = matureVecs.mean_vec;
std_vec = matureVecs.std_vec;
ages = matureVecs.ages;
inh_vec = matureVecs.inh_vec;
HVC_vec = matureVecs.HVC_vec;
HVCisi_m_vec = matureVecs.HVCisi_m_vec;
HVCisi_s_vec = matureVecs.HVCisi_s_vec;
HVCnum_m_vec = matureVecs.HVCnum_m_vec;
HVCnum_s_vec = matureVecs.HVCnum_s_vec;

%Cycle through each age to calculate inputs and simulate the neuron
for i = 1:numel(ages)
    
    %Use loop pointer to select out weight values from the HVC input vectors
    pruning = prune_vec(i);
    mean_w = mean_vec(i);
    var_w = std_vec(i)^2;
    
    %Generate lognormal mu and sigma from the HVC input vectors values
    mu = 2 * log(mean_w) - 0.5 * log(var_w + mean_w^2);
    sig = sqrt(log(var_w + mean_w^2) - 2 * log(mean_w));
    
    %HVC-RA synaptic weights from lognormal distribution constants
    %W_HVC = exp(randn(1,N_HVC) * sig + mu); %exponential distribution
    W_HVC = lognrnd(mu, sig, 1, N_HVC); %lognormal distribution
    
    %Reports neuron #, step # and corresponding synapse weight mean and std to the display
    %disp([pp, p, mean(W_HVC(W_HVC>0)), std(W_HVC(W_HVC>0))]);
    
    % Pruning
    W_HVC = W_HVC .* (rand(size(W_HVC)) >= pruning);
    
    %HVC-mediated inhibition
    %inh_HVC = mean(W_HVC) * inh_HVC_ratio;
    inh_HVC = mean(W_HVC) * inh_HVC_ratio * inh_vec(i); %Include the age-related change
    
    if staticHVC
        %Generate input from HVC that is exactly the same from rendition to
        %rendition
        
        %Template for one burst of HVC (5 spikes in 10ms by default; updated version scales HVC ISI by age
        %HVC_burst_template = repmat([1, zeros(1, round((2-dt)/dt))], 1 , 5);
        HVC_burst_template = repmat([1, zeros(1, round((HVC_vec(i)-dt)/dt))], 1 , 5);
        
        % Generating the whole HVC spikes (modulated by the corresponding rate)
        HVC_spikes = [];
        for j = 1:N_HVC
            HVC_spikes = [ HVC_spikes, HVC_burst_template * W_HVC(j)];
        end
        
        % Equalizing the length to t
        HVC_spikes = [HVC_spikes, zeros(size(t))];
        HVC_spikes = HVC_spikes(1: length(t));
        
        %Replicate the pattern for the full set of renditions
        HVC_spikes_array = repmat(HVC_spikes, numTrials, 1);

    else
        %Generate an HVC spiking pattern that has some pre-specified degree
        %of variability in the burst structure.
        
        % Constants/unpacking
%         m_numSpikes = 3.86; %from TMO estimate of Okubo et al 2015 data
%         s_numSpikes = 1.26; %from TMO estimate of Okubo et al 2015 data
%         m_ISI = 2.22; %from TMO estimate of Okubo et al 2015 data
%         s_ISI = 2.85; %from TMO estimate of Okubo et al 2015 data
        
        m_numSpikes = HVCnum_m_vec(i); %from TMO estimate of Okubo et al 2015 data
        s_numSpikes = HVCnum_s_vec(i); %from TMO estimate of Okubo et al 2015 data
        m_ISI = HVCisi_m_vec(i); %from TMO estimate of Okubo et al 2015 data
        s_ISI = HVCisi_s_vec(i); %from TMO estimate of Okubo et al 2015 data
        
        burstWidth = 10; % burst width (ms); 
        
        %Generate lognormal mu and sigma from the HVC spikes/burst
%         mu_numSp = 2 * log(m_numSpikes) - 0.5 * log(s_numSpikes + m_numSpikes^2);
%         sig_numSp = sqrt(log(s_numSpikes + m_numSpikes^2) - 2 * log(m_numSpikes));
        mu_numSp = m_numSpikes;
        sig_numSp = s_numSpikes;
        
        %Generate lognormal mu and sigma from the HVC ISI descriptive statistics
         mu_ISI = 2 * log(m_ISI) - 0.5 * log(s_ISI + m_ISI^2);
         sig_ISI = sqrt(log(s_ISI + m_ISI^2) - 2 * log(m_ISI));
%         mu_ISI = m_ISI;
%         sig_ISI = s_ISI;
        
        % Cycle through each rendition
        HVC_spikes_array = [];
%         numSpikesACC = [];
%         ISIACC = [];
%         clrIdx = randperm(N_HVC);
        for j = 1:numTrials
            
            %Save time by pre-selecting the number of spikes for the trial
            numSpikes = round(normrnd(mu_numSp, sig_numSp, [1, N_HVC]));        %normal distribution
            numSpikes(numSpikes<1) = 1;
            
            % Generate the HVC spikes for one trip through the synfire-chain (modulated by the corresponding rate)
            HVC_spikes = [];
            for k = 1:N_HVC
                
                %Processing time shortcut
                if W_HVC(k) == 0
                    %If the synapse weight is zero, just fill the template with zeros now
                    HVC_burst_template = zeros(1, burstWidth/dt);
                    
                else
                    %Retrive the correct number of spikes and ISIs as defined by distributions
                    %numSpikes = max([1, round(normrnd(mu_numSp, sig_numSp, 1))]);        %normal distribution
                    ISIs = lognrnd(mu_ISI, sig_ISI, 1, numSpikes(k));      %lognormal distribution
                    
                    %Assemble single burst
                    HVC_burst_template = [];
                    for l = 1:numSpikes(k)
                        timebins = round((ISIs(l)-dt)/dt);
                        spike = [1, zeros(1, timebins)];
                        
                        HVC_burst_template = [HVC_burst_template, spike];
                    end
                    
                    %Set burst width to defined value
                    HVC_burst_template = [HVC_burst_template, zeros(1, burstWidth/dt)];
                    HVC_burst_template = HVC_burst_template(1:burstWidth/dt);
                    
                    
                end
                
                %Assemble the spiking array from component bursts
                HVC_spikes = [HVC_spikes, HVC_burst_template * W_HVC(k)];
                %HVC_spikes = [HVC_spikes, HVC_burst_template * clrIdx(k)]; %Test to ID neuron boundaries
            end
            
            % Equalizing the length to t
            HVC_spikes = [HVC_spikes, zeros(size(t))];
            HVC_spikes = HVC_spikes(1:numel(t));
            
            %Copy current rendition to the output structure
            HVC_spikes_array(j,:) = HVC_spikes;
            
        end
    end
    
    %Visualize the spikes
%     figure(2001);
%     imagesc(HVC_spikes_array)
    
    %Simulate neuron given this input
%     [spike_times, Vm] = RA_spiking_neuron(HVC_spikes, inh_HVC, numTrials, t, numHVC, numLMAN);
    [spike_times, Vm] = RA_spiking_neuron(HVC_spikes_array, inh_HVC, numTrials, t, numHVC, numLMAN);
    
    %Print the rasters, if desired
    if raster && ((i==1) || ages(i)==62 || ages(i)==110 || (i == (length(ages))))
        
        figure(100);
        if i == 1
            %Plot the first 'day' of the simulation in subplot 1
            subplot(4, 1, 1);
            age = '42';
            
        elseif ages(i) == 62
            %Plot the first transition 'day' of the simulation in subplot 2
            subplot(4, 1, 2);
            age = '62';
            
        elseif ages(i) == 110
            %Plot the second transition 'day' of the simulation in subplot 3
            subplot(4, 1, 3)
            age = '110';
            
        elseif i == numel(ages)
            %Plot the last 'day' of the simulation in subplot 4
            subplot(4, 1, 4)
            age = '200';
        end
        
        %Streamline plotting code for rasters
        cla; hold on
        for j = 1:numel(spike_times)
            [tSx, tSy]  = rasterLine(spike_times{j}', (j - (lineSize/2)), (j + (lineSize/2)));
            line(tSx, tSy,'color','k');
        end
        xlim([0, t(end)]); ylim([0, j+1])
        ylabel(['Age: ' age 'd'])
        set(gca, 'Box', 'off', 'TickDir', 'out')
        
        set(gcf, 'Units', 'inches', 'Position', [7, 1, 7, 10])
    end
    
    %Calculate the spiketrain statistics for this timepoint
    %Firing rate
    FR(i) = makeFR(spike_times, t_stop);
    
    %Burstiness
    BF(i) = makeBurstFract(spike_times);
    
    %Sparseness
    sparse(i) = makeSparse(spike_times, t_stop, PSTH_bin_size);
    
    %New correlation metric
    spikeCorr(i) = makeCorr(spike_times, t_stop/1000);
    
end
drawnow

%Copy statistic arrays to the output structure
spikeStats.ages = ages;
spikeStats.FR = FR;
spikeStats.BF = BF;
spikeStats.sparse = sparse;
spikeStats.spikeCorr = spikeCorr;

%Dump the paramaters to the output file
params = [];


function FR = makeFR(RA_spiketimes, motifLength)
%This function handles the calculation of burst fraction for a given set of already aligned (or not) spikes.
%Plotting is handled by the calling function. Ready-to-print spike times for the desired cells are passed in
%handles.alignedSpikes.

%Initialize variables
FR = [];

%Calculate number of spikes across motifs
t = cell2mat(cellfun(@(x) length(x), RA_spiketimes, 'UniformOutput', 0));

%Calculate mean firing rate (in Hz) by dividing by recording time
FR = mean(t./(motifLength/1000));

function BF = makeBurstFract(RA_spiketimes)
%This function handles the calculation of burst fraction for a given set of already aligned (or not) spikes.
%Plotting is handled by the calling function. Ready-to-print spike times for the desired cells are passed in
%handles.alignedSpikes.

%Initialize variables
BF = [];

%Set burst constants
burstThresh = 150; %in Hz
thresh = 1/burstThresh; %convert to seconds

%Calculate per cell
t = cellfun(@(x) diff(x./1000), RA_spiketimes, 'UniformOutput', 0);

for j = 1:numel(t)
    %Get indices of spikes in burst
    HF_index = find(t{j}<thresh);
    spikes_adj = length(find(diff(HF_index)==1));
    
    if ~isempty(spikes_adj) || spikes_adj ~= 0
        %Calculate the total number of spikes in bursts
        HF_spikes = 2*length(HF_index)-spikes_adj;
        
        %Calculate the burst fraction for the motif
        BF(j) = HF_spikes/(length(RA_spiketimes{j}));
    else
        BF(j) = 0;
        
    end
    
end
BF = nanmean(BF);

function sparse = makeSparse(RA_spiketimes, motifLength, binSize)
%This function handles the calculation of sparseness for a given set of already aligned (or not) spikes.
%Plotting is handled by the calling function. Ready-to-print spike times for the desired cells are passed in
%handles.alignedSpikes.

%Initialize variables
binCounts = [];
binsFR = [];

%Set PSTH bins
edges = 0:binSize:motifLength;

%Pool spike times (for each cell) across motifd
sumSpikeTimes = cell2mat(RA_spiketimes);

%Bin spikes
[binCounts,~] = histcounts(sumSpikeTimes,edges);

%Convert units to mean firing rate
numCellRend = numel(RA_spiketimes);
binsFR = binCounts./(numCellRend*binSize/1000);

%Convert to Spiking PDF by normalizing
pdfFR = binsFR./sum(binsFR);
L = log(pdfFR); L(isinf(L)) = 0; %Correct for log(0)
sparse = 1 + sum(pdfFR.*L)/log(length(pdfFR));

if isnan(sparse)
    sparse = 0;
end

function spikeCorr = makeCorr(RA_spiketimes, motifLength)
%This function handles the calculation of precision/correlation for a given set of already aligned (or not) spikes.
%Plotting is handled by the calling function. Ready-to-print spike times for the desired cells are passed in
%handles.alignedSpikes.

%Initialize variables
spikeCorr = [];
fs = 44150;

%Set correlation constants
gaussWidth=0.008; %in seconds
binSize = 1/fs;
trainLength = floor((motifLength + eps)/ binSize + 1);

%Generate gaussian
sigma = gaussWidth / sqrt(2);
x = [gaussWidth*-4:1/fs:gaussWidth*4];
gauss = (1/sqrt(2*pi)*sigma)*exp(-x.^2/(2*sigma^2));

%Convert spike times to a matrix of spike train signals (cov w/ gaussian)
sigs = cellfun(@(x) makeTrainSigs(x, gauss, trainLength, binSize), RA_spiketimes, 'UniformOutput', 0);
sigMat = cell2mat(sigs);

%Calculate the Correlation matrix
corrMat = corrcoef(sigMat);

%Report only the mean correlation (minus the autocorr, i.e. the diag)
sum_corr = nansum(nansum(corrMat))-nansum(diag(corrMat));
[m,n] = size(corrMat);
spikeCorr = sum_corr/((m*n)-m);

function sig = makeTrainSigs(x, gauss, trainLength, binSize)
%This function converts an array of spiketimes (relative to a motif start) and converts it into an analog signal. Mainly use
%for calulating pairwise correlations

%Create a binary spike train
spikeNdx = floor((x/1000 + eps) / binSize + 1);
binSpikes = zeros(trainLength, 1);
binSpikes(spikeNdx) = 1;

%Convolve binary train with the passed gaussian
sig  = conv(binSpikes, gauss);
sig = sig(floor(length(gauss)/2):end - floor((length(gauss)+1)/2)); %trim convolution artifact


