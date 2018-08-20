function simpleRAmodel
% function simpleRAmodel
% written by TMO 2017; last modified 10/09/17
%
% This script simulates a model RA neuron (i.e., it models changes in spiking as a
% function of HVC, LMAN, and local inhibition. The code is based on the
% "RA_toy_model_lognormal_orig" script used in Garst-Orozco et al (2015)
% and incorporates some lessons learned along the way.
%
% As a first pass, this will be a stripped down version of the
% matureRA_tvarHVC.m model.
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
numTrials = 10; %Number of trials to run per age per neuron
N_HVC = 100; %Number of HVC neurons in complete pool
numLMAN = 2; %LMAN input multiplier
lineSize = 0.8; %raster plot line height
PSTH_bin_size = 20; %bin size in (ms)
inh_HVC_ratio = 1600;
numHVC = 1;

%Simulation time vector
t_stop = 1000;  % length of each trial (ms);
dt = 0.2; %time step (ms)
t = dt:dt:t_stop;
nodeSpace = 10; %spacing of HVC nodes in the model

% %Setup the vectors defining the maturation of connectivity
% prune_vec = matureVecs.prune_vec;
% mean_vec = matureVecs.mean_vec;
% std_vec = matureVecs.std_vec;
% ages = matureVecs.ages;
% inh_vec = matureVecs.inh_vec;
% HVC_vec = matureVecs.HVC_vec;
% HVCisi_m_vec = matureVecs.HVCisi_m_vec;
% HVCisi_s_vec = matureVecs.HVCisi_s_vec;
% HVCnum_m_vec = matureVecs.HVCnum_m_vec;
% HVCnum_s_vec = matureVecs.HVCnum_s_vec;

%Test age point... replace with an array when ready
ages = 70;
paradigm = 'tutored';

%Cycle through each age to calculate inputs and simulate the neuron
for i = 1:numel(ages)
    %Retrieve the parameters for generating an HVC pool at this age
    [mu_ISI, sig_ISI, mu_numBurst, sig_numBurst] = genHVCparams(ages(i));
    

    
    %Retrieve age-dependent synapse parameters (i.e., strength and number)
    [mu_SF, sig_SF, prune] = genRAparams(ages(i), paradigm);
    
    %HVC-RA synaptic weights from lognormal distribution constants
    %W_HVC = exp(randn(1,N_HVC) * sig + mu); %exponential distribution
    W_HVC = lognrnd(mu_SF, sig_SF, 1, N_HVC); %lognormal distribution
    
    %Reports neuron #, step # and corresponding synapse weight mean and std to the display
    %disp([pp, p, mean(W_HVC(W_HVC>0)), std(W_HVC(W_HVC>0))]);
    
    % Random pruning (i.e., randomly zero the weights)
    W_HVC = W_HVC .* (rand(size(W_HVC)) >= prune);
    
    %Generate index of non-zero weights (for later use)
    wIdx = find(W_HVC ~= 0);
    
    %HVC-mediated inhibition
    inh_HVC = mean(W_HVC) * inh_HVC_ratio;
    %inh_HVC = mean(W_HVC) * inh_HVC_ratio * inh_vec(i); %Include the age-related change
    
    %Generate a new set of spike trans for each trial rendition (but keep
    %the weights (and thus pruning) consistent between rensitions)
    for j = 1:numTrials
        %Generate the HVC pool from the array
        spikeMat = genHVCpool(N_HVC, mu_ISI, sig_ISI, mu_numBurst, sig_numBurst, dt, nodeSpace);
        
        %Generate the weighted spike times
        for k = wIdx
            spikeMat(k,:) = W_HVC(k) .* spikeMat(k,:);
        end
        
        %Sum these spike weights across neurons for net input at ea h moment of time
        HVC_spikes(j,:) =  sum(spikeMat,1);
        
    end
end


%Visualize the spikes
%     figure(2001);
%     imagesc(HVC_spikes_array)

%Simulate neuron given this input
%     [spike_times, Vm] = RA_spiking_neuron(HVC_spikes, inh_HVC, numTrials, t, numHVC, numLMAN);
[spike_times, Vm] = RA_spiking_neuron(HVC_spikes, inh_HVC, numTrials, t, numHVC, numLMAN);

%Print the rasters, if desired
% if raster && ((i==1) || ages(i)==62 || ages(i)==110 || (i == (length(ages))))
if raster   
    figure(100);
%     if i == 1
%         %Plot the first 'day' of the simulation in subplot 1
%         subplot(4, 1, 1);
%         age = '42';
%         
%     elseif ages(i) == 62
%         %Plot the first transition 'day' of the simulation in subplot 2
%         subplot(4, 1, 2);
%         age = '62';
%         
%     elseif ages(i) == 110
%         %Plot the second transition 'day' of the simulation in subplot 3
%         subplot(4, 1, 3)
%         age = '110';
%         
%     elseif i == numel(ages)
%         %Plot the last 'day' of the simulation in subplot 4
%         subplot(4, 1, 4)
%         age = '200';
%     end
    
    %Streamline plotting code for rasters
    cla; hold on
    for j = 1:numel(spike_times)
        [tSx, tSy]  = rasterLine(spike_times{j}', (j - (lineSize/2)), (j + (lineSize/2)));
        line(tSx, tSy,'color','k');
    end
    xlim([0, t(end)]); ylim([0, j+1])
%     ylabel(['Age: ' age 'd'])
    set(gca, 'Box', 'off', 'TickDir', 'out')
    
%     set(gcf, 'Units', 'inches', 'Position', [7, 1, 7, 10])
end

%Calculate the spike train stats
if false
    %Calculate the spiketrain statistics for this timepoint
    %Firing rate
    FR(i) = makeFR(spike_times, t_stop);
    
    %Burstiness
    BF(i) = makeBurstFract(spike_times);
    
    %Sparseness
    sparse(i) = makeSparse(spike_times, t_stop, PSTH_bin_size);
    
    %New correlation metric
    spikeCorr(i) = makeCorr(spike_times, t_stop/1000);
    
    drawnow
    
    %Copy statistic arrays to the output structure
    spikeStats.ages = ages;
    spikeStats.FR = FR;
    spikeStats.BF = BF;
    spikeStats.sparse = sparse;
    spikeStats.spikeCorr = spikeCorr;
    
    %Dump the paramaters to the output file
    params = [];
end

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










