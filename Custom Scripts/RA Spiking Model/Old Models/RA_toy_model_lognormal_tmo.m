function RA_toy_model_lognormal_tmo
% RA_toy_model_lognormal_tmo.m
% edited by TMO 05-23-2016
%
% This script is a commented and edited version of Baktash's script
% RA_toy_model_lognormal_orig.m that Bence provided sometime in 2016.

% clear
raster = 1; %boolean flag for plotting output

%%%%%%%%%%%%%%%%%%%%%%%%%
%Set simulation constants
%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the number of RA neurons to simulate
pool = 100;

%Define the number of trials for each neuron
trials = 100; 

%Define the length of the simulated "song" in a trial
t_stop = 1000;  % length of each trial (ms);
dt = 0.2; %time step (ms)
t = dt : dt : t_stop; %time array

%The section defines the developmental changes in HVC input to RA. 'mean'
%and 'std' define the synapse weight distributions (in pA?); 'prune' is the
%fraction of synapses pruned away. The values are based on the in vitro
%measurements, though interpolation and extrapolation is used to increase
%resolution. unknown why this many steps were chosen.
%
%Generate the pruning and weights vectors for the original simulation
%[prune_vec, mean_vec, std_vec] = originalVects;
%inh_vec = ones(size(prune_vec));

%Generate the pruning and weights vectors for the new simulation
[prune_vec, mean_vec, std_vec, ages, inh_vec, HVC_vec] = expVects;

%Force constant vectors
a = 6;
if a == 1
    %Constant pruning
    const = prune_vec(1);
    new_vec = const.*ones(size(prune_vec));
    prune_vec = new_vec;
    
elseif a == 2
    %Constant weights
    const = mean_vec(1);
    new_vec = const.*ones(size(mean_vec));
    mean_vec = new_vec;
    
    const = std_vec(1);
    new_vec = const.*ones(size(std_vec));
    std_vec = new_vec;
    
elseif a == 3
    %Constant inhibition
    const = inh_vec(1);
    new_vec = const.*ones(size(inh_vec));
    inh_vec = new_vec;
    
elseif a == 4
    %Constant HVC ISI
    const = HVC_vec(end); %Default to adult
    new_vec = const.*ones(size(HVC_vec));
    HVC_vec = new_vec;
end

%RA neuron parameters (time constants, voltages, input weights)
tau_s = 5;      % RA synapse time constant (ms)
tau_m = 20;     % RA membrane time constant (ms)
V_th = 20;      % RA firing threshold (mV)
V_r = 0;        % RA resting membrane voltage (everything is offset by +70mV)
t_ref = 1;      % Refractory period of RA (ms)
tau_NMDA = 100; % NMDA receptor time constant (ms)
c_Mg = 0.7;     % Mg ion concentration (mM)
LMAN_res = 1 * 130; % LMAN input resistance (in MOhms)
HVC_res = 4 * 130;  % HVC input resistance (in MOhms)

%Inhibition
% inh_HVC_ratio = 800;
inh_HVC_ratio = 400;
inh_LMAN = 0;

%LMAN input
LMAN_rate = 80; % Rate of LMAN input (Hz)

%Alternative LMAN inputs characteristics
noise_std = 0;
W_LMAN =  150/1000;   % LMAN-RA synaptic strength (pA) For NMDA
N_ratio = 0.9; % NMDA/(NMDA+AMP) ratio

% % noise_std = 0;
% % W_LMAN = 90;   % LMAN-RA synaptic strength (mV) For NMDA
% % inh = 150;
% % N_ratio = 0.85; % NMDA/(NMDA+AMP) ratio
% 
% % noise_std = 0;
% % W_LMAN = 220;   % LMAN-RA synaptic strength (mV) for AMPA
% % inh = 45;
% % N_ratio = 0; % NMDA/(NMDA+AMP) ratio(
% 
% % noise_std = 150;
% % inh = 100;
% % W_LMAN = 0;   % LMAN-RA synaptic strength (mV) for AMPA
% % N_ratio = 0; % NMDA/(NMDA+AMP) ratio

%HVC network size
N_HVC = 100; % Number of HVC neurons

%Matrix of HVC inputs to RA
W_HVC_mat = zeros(length(mean_vec),pool*N_HVC); 

% %Template for one burst of HVC (5 spikes in 10ms by default)
% HVC_burst_template = repmat([1, zeros(1, round((2-dt)/dt))], 1 , 5);

%
PSTH_bin_size = 20;

%Predefine the output structures
CC_vec =  zeros(size(prune_vec));
CC_vec_2 = zeros(size(prune_vec));
rate_vec = zeros(size(prune_vec));
rec = zeros(size(t)); %This looks like garbage

FR_vec = zeros(size(prune_vec));
BF_vec = zeros(size(prune_vec));
sparse_vec = zeros(size(prune_vec));
spikeCorr_vec = zeros(size(prune_vec));
I_vec = zeros(size(prune_vec));
HVCI_vec = zeros(size(prune_vec));
MAX_vec = zeros(size(prune_vec));

%%%%%%%%%%%%%%%%%%%%%%%%%
%Run simulation
%%%%%%%%%%%%%%%%%%%%%%%%%

%Cycle through the pool of RA neurons
for pp = 1:pool
    
    %Clear the output figures if rasters outputs are selected
    if (raster == 1)
        figure(2);
        clf
    end
    
    %Cycle through each of the pruning steps
    for p = 1:length(prune_vec)

        %Use loop pointer to select out weight values from the HVC input vectors
        pruning = prune_vec(p);
        mean_w = mean_vec(p);
        var_w = std_vec(p)^2;
        
        %Generate lognormal mu and sigma from the HVC input vectors values
        mu = 2 * log(mean_w) - 0.5 * log(var_w + mean_w^2);
        sig = sqrt(log(var_w + mean_w^2) - 2 * log(mean_w));
        
        %Seems to be unused?
        %rec_c = zeros(1,length(t) * trials);
        
        %HVC-RA synaptic weights from lognormal distribution constants
        %W_HVC = exp(randn(1,N_HVC) * sig + mu); %exponential distribution
        W_HVC = lognrnd(mu, sig,1,N_HVC); %lognormal distribution
        
        %Reports neuron #, step # and corresponding synapse weight mean and std to the display
        disp([pp, p, mean(W_HVC(W_HVC>0)), std(W_HVC(W_HVC>0))]);
        
        % Pruning
        W_HVC = W_HVC .* (rand(size(W_HVC)) >= pruning);
        
        %HVC-mediated inhibition
        %inh_HVC = mean(W_HVC) * inh_HVC_ratio;
        inh_HVC = mean(W_HVC) * inh_HVC_ratio * inh_vec(p); %Include the age-related change
        
        inh = inh_LMAN + inh_HVC;
        
        %Template for one burst of HVC (5 spikes in 10ms by default;
        %updated version scales HVC ISI by age
%       HVC_burst_template = repmat([1, zeros(1, round((2-dt)/dt))], 1 , 5);
        HVC_burst_template = repmat([1, zeros(1, round((HVC_vec(p)-dt)/dt))], 1 , 5);

        % Generating the whole HVC spikes (modulated by the corresponding rate)
        HVC_spikes = [];
        for i = 1 : N_HVC
            HVC_spikes = [ HVC_spikes, HVC_burst_template * W_HVC(i)];
        end
        
        % Equalizing the length to t
        HVC_spikes = [HVC_spikes, zeros(size(t))];
        HVC_spikes = HVC_spikes(1: length(t));
        RA_spikes = zeros(trials, length(t));
        
        %This is old code from the original, but not sure of it's
        %relevance... 
%         HVC_spikes = zeros(size(HVC_spikes));
%         HVC_spikes(10) = 0;
%         W_HVC = ones(size(W_HVC));
%         W_LMAN = 1;   % LMAN-RA synaptic strength (mV) For NMDA
%         inh = 0;
%         N_ratio = 1;
        
        prev_spike_time = 0;
        ISI = zeros(trials, length(t));
        recSum = zeros(trials,1);
        %Loop on trials from this particular neuron
        for i = 1 : trials
            
            % LMAN input to RA neuron (rate model of random input)
            LMAN_spikes = rand(size(t))< (LMAN_rate * dt / 1000);

            if i == 1
                LMAN_AMPA_input = zeros(size(t));
                LMAN_NMDA_input = zeros(size(t));
                HVC_input = zeros(size(t));
                V_RA = zeros(size(t));
            else
                LMAN_AMPA_input(1) = LMAN_AMPA_input(end);
                LMAN_NMDA_input(1) = LMAN_NMDA_input(end);
                HVC_input(1) = 1.5 * HVC_input(tt);
                V_RA(1) = V_RA(end);
            end
            
            %Loop on time through the simulation
            ref = 0;
            prev_tt = 1;
            for tt = 2 : length(t)
                
                %Mg-block
                x_Mg = 1/(1 + exp(- 0.062 * (V_RA(tt-1)-70)) * c_Mg/3.57);
                
                %LMAN input rate of change
                dLMAN_AMPA_input_dt = (- LMAN_AMPA_input(tt - 1) + tau_s * (1-N_ratio) * W_LMAN * LMAN_spikes(tt)/dt)/tau_s;
                dLMAN_NMDA_input_dt = (- LMAN_NMDA_input(tt - 1) + tau_NMDA * N_ratio * x_Mg * W_LMAN * LMAN_spikes(tt)/dt)/tau_NMDA;
                
                %HVC input rate of change
                dHVC_input_dt = ( - HVC_input(tt - 1) + tau_s * HVC_spikes(tt)/dt)/tau_s;
                
                %LMAN input
                LMAN_AMPA_input(tt) = LMAN_AMPA_input(tt - 1) + dt * dLMAN_AMPA_input_dt;
                LMAN_NMDA_input(tt) = LMAN_NMDA_input(tt - 1) + dt * dLMAN_NMDA_input_dt;
                
                %HVC input
                HVC_input(tt) = HVC_input(tt - 1) + dt * dHVC_input_dt;
                
                %Randomized noise added to the simulation
                I_noise = noise_std/sqrt(dt) * randn;
                
                %RA neuron voltage rate of change summations for all inputs
                dV_RA_dt = (-V_RA(tt-1) + LMAN_AMPA_input(tt) * LMAN_res + LMAN_NMDA_input(tt) * LMAN_res + HVC_input(tt) * HVC_res - inh + I_noise)/tau_m;
                
                if ref<=0
                    V_RA(tt) = V_RA(tt-1) + dt * dV_RA_dt;
                else
                    V_RA(tt) = 0;
                end
                ref = ref - dt;
                
                %Spiking conditions for RA neurons
                if (V_RA(tt) >=V_th) && (ref<=0)
                    RA_spikes(i,tt) = 1;
                    V_RA(tt) = V_r;
                    ref = t_ref;
                    ISI(i, prev_tt : tt) = t(tt) - prev_spike_time;
                    prev_spike_time = t(tt);
                    prev_tt = tt;
                end
                
            %This appears to do nothing important...???    
            rec(tt) = HVC_input(tt) * HVC_res;
            
            end
            
            %Total HVC excitation (during this particular trial)
            recSum(i) = mean(rec);
            
%             figure(100); clf
%             plot(1:tt, V_RA, 'r')
%             a = 1;
            
            %Collect spike times for later stats
%             RAst = RA_spikes(i,:) .* t;
%             RAst(RAst == 0) = [];
%             RA_spike_times_perm{i} = RAst;
            RA_spike_times = RA_spikes(i,:) .* t;
            RA_spike_times(RA_spike_times == 0) = [];
            RA_spike_times_perm{i} = RA_spike_times;

            % RA Raster Plot (only plot the first 100 trials)
            if (i<=100) && (raster == 1) && ((p==1) || ages(p)==62 || ages(p)==110 || (p == (length(prune_vec))))
%                 RA_spike_times = RA_spikes(i,:) .* t;
%                 RA_spike_times(RA_spike_times == 0) = [];
%                 RA_spike_times_perm{i} = RA_spike_times;
                
                 %Raster plot figure switch
%                 if p == 1
%                     %Plot the first 'day' of the simulation in Figure 2
%                     figure(2);
%                     xlabel(num2str(pruning));
%                 else
%                     %Plot the all subsequent 'days' in Figure 3
%                     figure(3);
%                     xlabel(num2str(pruning));
%                 end;
                
                figure(2);
                if p == 1
                    %Plot the first 'day' of the simulation in subplot 1
                    subplot(4, 1, 1)
                    age = '42';
                    
                elseif ages(p) == 62
                    %Plot the first transition 'day' of the simulation in subplot 2
                    subplot(4, 1, 2)
                    age = '62';
                    
                elseif ages(p) == 110
                    %Plot the second transition 'day' of the simulation in subplot 3
                    subplot(4, 1, 3)
                    age = '110';
                    
                elseif p == length(prune_vec)
                    %Plot the last 'day' of the simulation in subplot 4
                    subplot(4, 1, 4)
                    age = '200';
                end

                %Streamline plotting code (Added 5/30/16)
                hold on
                [tSx, tSy]  = rasterLine(RA_spike_times', (i-1)*3, (i-1)*3+2);
                line(tSx, tSy,'color','k');
                ylabel(['Age: ' age 'd'])
                
                set(gca, 'Box', 'off', 'TickDir', 'out')
                set(gcf, 'Units', 'inches', 'Position', [7, 1, 7, 10])
                
            end
            
        end
        
        drawnow
        
%         %Calulcate PSTH
%         bin_num = round(PSTH_bin_size/dt);
%         for i = 1 : trials
%             RA_spikes_bined(i,:) = sum(reshape(RA_spikes(i,:), bin_num, round(length(t)/bin_num)),1);
%         end;
%         
%         %Calculate firing rate for each trial from PSTH
%         RA_rate_mean = 1000* mean(RA_spikes_bined,1)/PSTH_bin_size;
%         RA_rate_std =  std(1000* RA_spikes_bined /(PSTH_bin_size),1)/sqrt(PSTH_bin_size);
%         
%         %Smoothing for each trial (Note that this is different from the method above)
%         RA_rate_smooth = 1000./ISI;
%         RA_rate_smooth(isnan(RA_rate_smooth)) = 0;
%         RA_rate_smooth(isinf(RA_rate_smooth)) = 0;
%         
%         %Old way that was commented out
%         % RA_rate_mean = 1000 * reshape(sum(RA_spikes,1), bin_num, round(length(t)/bin_num))/ (trials * dt);
%         % RA_rate_mean = mean(RA_rate_mean,1);
%         %      = smoothts(1000./ISI,'g',1000,10/dt);
%         
%         %Correlation coefficients based on the ISI
%         C = corrcoef(RA_rate_smooth');
%         %C = corrcoef(RA_spikes_bined');
%         
%         C = C .* (1 - eye(trials));
%         C(isnan(C)) = 0;
%         CC = sum(C(:)) / (trials^2 - trials);
% %         if RA_rate_mean <20
% %             CC = CC_vec(p);
% %         end;
% 
%         %Calculate the running average for the correlations
%         CC_vec(p) =  ((pp - 1) * CC_vec(p) + CC)/pp;
% %         CC_vec_2(p)= ((pp - 1) * CC_vec_2(p) + CC^2)/pp;
%         
%         %Calculate the runnning average for the firing rates
%         rate_vec(p) = ((pp - 1) * rate_vec(p) + mean(RA_rate_mean))/pp;
        
        %Seems like this can output the distribution of HVC-RA synapse weights
%         [mean(LMAN_AMPA_input), std(LMAN_AMPA_input), 1000*sum(sum(RA_spikes))/t_stop/trials]
%         W_HVC_mat(p,( (pp-1) * N_HVC + 1): pp * N_HVC) = W_HVC;
%         W = W_HVC_mat(p,:);
%         figure(4);
%         hand = subplot(length(prune_vec),1,p);
%         hist(W(W>0),20);
%         set(hand,'xlim', [0, 0.3]);
%         disp([p, mean(W(W>0)), std(W(W>0))]);

        %Firing rate
        FR(p) = makeFR(RA_spike_times_perm, t_stop);
        
        %Burstiness
        BF(p) = makeBurstFract(RA_spike_times_perm);

        %Sparseness
        sparse(p) = makeSparse(RA_spike_times_perm, t_stop, PSTH_bin_size);

        %New correlation metric
        spikeCorr(p) = makeCorr(RA_spike_times_perm, t_stop/1000);
        
        %Inhibition
        I(p) = inh;
        
        %HVC ISI
        HVCI(p) = mean(recSum);
        
        %Sum of the weights (should be proportional to the MAX current from JGO 2015)
        MAX(p) = sum(W_HVC(W_HVC>0));
        
        %Running vectors
        %Calculate the running average for the correlations
        FR_vec(p) =  ((pp - 1) * FR_vec(p) + FR(p))/pp;
        BF_vec(p) =  ((pp - 1) * BF_vec(p) + BF(p))/pp;
        sparse_vec(p) =  ((pp - 1) * sparse_vec(p) + sparse(p))/pp;
        spikeCorr_vec(p) =  ((pp - 1) * spikeCorr_vec(p) + spikeCorr(p))/pp;
        I_vec(p) =  ((pp - 1) * I_vec(p) + I(p))/pp;
        HVCI_vec(p) =  ((pp - 1) * HVCI_vec(p) + HVCI(p))/pp;
        MAX_vec(p) =  ((pp - 1) * MAX_vec(p) + MAX(p))/pp;
    end
    

% %Plot the running statistics for the simulation
% figure(1); clf
% 
% %Plot the pruning fraction vs mean firing rate
% plot(prune_vec, rate_vec,'r'); hold on;
% 
% %Plot the pruning fraction vs cross-correlation
% plot(prune_vec, CC_vec*100,'b');
% drawnow
% 
% xlabel('Pruning Fraction')
% ylabel('FR or CC')
% title('Simulation Statistics')

%Plot all
figure(57); clf

subplot(7,1,1)
[hAx, hLine1, hLine2] = plotyy(ages, prune_vec, ages, mean_vec);
ylabel(hAx(1), 'Pruning Fraction')
ylabel(hAx(2), 'Synapse Weight')
title('Simulation Statistics')
xlim(hAx(1), [35, 210])%; ylim(hAx(1), [0, 1])
xlim(hAx(2), [35, 210])
set(hAx(1), 'Box', 'off', 'TickDir', 'out')
set(hAx(2), 'Box', 'off', 'TickDir', 'out')

subplot(7,1,2)
[hAx, hLine1, hLine2] = plotyy(ages, inh_vec, ages, I_vec);
ylabel(hAx(1), 'Inhibition Scale')
ylabel(hAx(2), 'Inhibition Sum')
xlim(hAx(1), [35, 210])%; ylim(hAx(1), [0, 1])
xlim(hAx(2), [35, 210])
set(hAx(1), 'Box', 'off', 'TickDir', 'out')
set(hAx(2), 'Box', 'off', 'TickDir', 'out')

subplot(7,1,3)
[hAx, hLine1, hLine2] = plotyy(ages, HVC_vec, ages, HVCI_vec);
ylabel(hAx(1), 'HVC ISI (ms)')
ylabel(hAx(2), 'HVC Sum')
xlim(hAx(1), [35, 210])%; ylim(hAx(1), [0, 1])
xlim(hAx(2), [35, 210])
set(hAx(1), 'Box', 'off', 'TickDir', 'out')
set(hAx(2), 'Box', 'off', 'TickDir', 'out')

subplot(7,1,4)
plot(ages, FR_vec,'b')
ylabel('Firing Rate (Hz)')
xlim([35, 210])
set(gca, 'Box', 'off', 'TickDir', 'out')

subplot(7,1,5)
plot(ages, BF_vec,'r')
ylabel('Burst Fraction')
xlim([35, 210]); ylim([0,1])
set(gca, 'Box', 'off', 'TickDir', 'out')

subplot(7,1,6)
plot(ages, sparse_vec,'g')
ylabel('Sparseness')
xlim([35, 210])
set(gca, 'Box', 'off', 'TickDir', 'out')

subplot(7,1,7)
plot(ages, spikeCorr_vec,'k');
xlabel('Ages (d)')
ylabel('Correlation')
xlim([35, 210]); ylim([0,1])
set(gca, 'Box', 'off', 'TickDir', 'out')

set(gcf, 'Units', 'inches', 'Position', [1, 1, 5, 12])
drawnow

% figure(59)
% plot(ages, MAX_vec,'k');
% xlabel('Ages (d)')
% ylabel('MAX')
% xlim([35, 210]);
% set(gca, 'Box', 'off', 'TickDir', 'out')
% drawnow

% %Plot ages vs stats
% figure(5); clf
% 
% subplot(3,1,1)
% [hAx, hLine1, hLine2] = plotyy(ages, prune_vec, ages, mean_vec);
% ylabel(hAx(1), 'Pruning Fraction')
% ylabel(hAx(2), 'Synapse Weight')
% 
% subplot(3,1,2)
% plot(ages, rate_vec,'r'); drawnow
% ylabel('FR')
% title('Simulation Statistics')
% 
% subplot(3,1,3)
% plot(ages, CC_vec*100,'b'); drawnow
% xlabel('Ages (d)')
% ylabel('CC')

end

function [prune_vec, mean_vec, std_vec, ages, inh_vec, HVC_vec] = expVects
% Empirical data from JGO (TMO calculated on 5/31 -- these are the descriptive stats)
%           T40       I40        T60       I60      T90       I90      T200      I200
%sfMean = 29.5219   29.9044   49.6013   45.3759   73.6766   72.9715   79.1057   66.9595
%sfStd  = 21.7046   20.1472   30.0967   30.4014   60.3826   72.7731   43.5912   57.0972
%inEst  = 18.6754   18.1951   26.3587   30.8637   10.8243   37.1721   10.0262   23.1487

%Steps in segments
% s1 = 11; %2d spacing
% s2 = 25; %2d spacing
% s3 = 19; %5d spacing

s1 = 6; %4d spacing
s2 = 8; %4d spacing
s3 = 12; %10d spacing

%Linear scaling function for the pruning
%The basic premise is that the value of rho should be proportional to the
%extent of synapse pruning such that
%
%min inputs => rho==1   all synapses pruned
%max inputs => rho==0   no synapses pruned
% ((in_max-in_min)-# ins)/in_max
in_max = 50;
in_min = 0;

%Sparse datapoints for d60, d110, and d200
sp_ages = [42, 62, 90, 200];

%Tutored
sp_Mean = [30, 50, 75, 80]./1000;
sp_Std  = [22, 30, 60, 44]./1000;
%sp_rawPrun = [18, 26.4, 10.8, 10];
sp_HVC = [6, 4, 2, 2]; %age-dependent decrease in HVC ISI (ms)
%sp_HVC = [2, 2, 2, 2]; %age-dependent decrease in HVC ISI (ms)
%sp_inh = [1, 1, 0.625, 0.625]; % age-dependent inhibition based on Sakaguchi 1996
sp_inh = [1, 2, 2, 2]; % age-dependent inhibition based on Olveczky 2011

%Isolate
% sp_Mean = [29.9, 45.4, 73, 67]./1000;
% sp_Std  = [20.1, 30.4, 72.8, 57.1]./1000;
sp_rawPrun = [18, 31, 37, 23];
% sp_inh = [1, 1, 0.625, 0.625]; % age-dependent inhibition based on Sakaguchi 1996

sp_Prun = ((in_max-in_min)-sp_rawPrun)./in_max; %Linear scaling

%Interpolate between points
r1 = linspace(sp_ages(1), sp_ages(2), s1); r2 = linspace(sp_ages(2), sp_ages(3), s2); r3 = linspace(sp_ages(3), sp_ages(4), s3);
ages = [r1, r2(2:end), r3(2:end)];

r1 = linspace(sp_Mean(1), sp_Mean(2), s1); r2 = linspace(sp_Mean(2), sp_Mean(3), s2); r3 = linspace(sp_Mean(3), sp_Mean(4), s3);
mean_vec = [r1, r2(2:end), r3(2:end)];

r1 = linspace(sp_Std(1), sp_Std(2), s1); r2 = linspace(sp_Std(2), sp_Std(3), s2); r3 = linspace(sp_Std(3), sp_Std(4), s3);
std_vec = [r1, r2(2:end), r3(2:end)];

r1 = linspace(sp_Prun(1), sp_Prun(2), s1); r2 = linspace(sp_Prun(2), sp_Prun(3), s2); r3 = linspace(sp_Prun(3), sp_Prun(4), s3);
prune_vec = [r1, r2(2:end), r3(2:end)];

r1 = linspace(sp_inh(1), sp_inh(2), s1); r2 = linspace(sp_inh(2), sp_inh(3), s2); r3 = linspace(sp_inh(3), sp_inh(4), s3);
inh_vec = [r1, r2(2:end), r3(2:end)];

r1 = linspace(sp_HVC(1), sp_HVC(2), s1); r2 = linspace(sp_HVC(2), sp_HVC(3), s2); r3 = linspace(sp_HVC(3), sp_HVC(4), s3);
HVC_vec = [r1, r2(2:end), r3(2:end)];

function [prune_vec, mean_vec, std_vec] = originalVects
%Original vectors from the JGO paper

%weights
max_mean = 75/1000;
min_mean = 50/1000;
max_std = 80 /1000;
min_std = 30/1000;

%pruning fractions
min_prune = 0.1; %i.e., 90% of possible inputs intact
max_prune = 1 - (1 - min_prune) * (10/24);

steps = 15;

mean_steps = (max_mean - min_mean)/steps;
std_steps = (max_std - min_std)/steps;
prune_steps = (max_prune - min_prune)/steps;

mean_vec = [min_mean : mean_steps : max_mean];
std_vec = [min_std : std_steps : max_std];
prune_vec = [min_prune : prune_steps : max_prune];

prune_vec_head = [0 : prune_steps : min_prune];
prune_vec_tail = [max_prune + prune_steps : prune_steps : 0.8];
head_length = length(prune_vec_head);
tail_length = length(prune_vec_tail);

mean_vec_head = [min_mean - mean_steps * head_length : mean_steps : min_mean - mean_steps];
std_vec_head = [min_std - std_steps * head_length : std_steps : min_std - std_steps];
mean_vec_tail = [max_mean + mean_steps : mean_steps : max_mean + mean_steps * tail_length];
std_vec_tail = [max_std + std_steps : std_steps : max_std + std_steps * tail_length];

%these are the fully assembled HVC input vectors. as written, these produce
%24-element vectors (approx, ~2 days/element) though I see no justification for that choice.
prune_vec = [prune_vec_head , prune_vec, prune_vec_tail];
mean_vec = [mean_vec_head, mean_vec, mean_vec_tail];
std_vec = [std_vec_head, std_vec, std_vec_tail];

function FR = makeFR(RA_spiketimes, motifLength)
%This function handles the calculation of burst fraction for a given set of already aligned (or not) spikes.
%Plotting is handled by the calling function. Ready-to-print spike times for the desired cells are passed in
%handles.alignedSpikes. 

%Initialize variables
FR = [];

%Determine which cells are selected for analysis
% indx = get(handles.list_cellSelect,'Value');
% selectCells = handles.cellList(indx);

%Calculate and plot PSTH per cell
% cellId = getFieldVectorCell(handles.filtSync, 'cellname');
% for i = 1:numel(selectCells)
    %Generate per-cell logical mask
%     mask = strcmp(selectCells{i}, cellId);
    
    %Calculate number of spikes across motifs
%     s = handles.alignedSpikes(mask);
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

%Determine which cells are selected for analysis
% indx = get(handles.list_cellSelect,'Value');
% selectCells = handles.cellList(indx);

%Calculate per cell
% cellId = getFieldVectorCell(handles.filtSync, 'cellname');
% for i = 1:numel(selectCells)
    %Generate per-cell logical mask
%     mask = strcmp(selectCells{i}, cellId);
    
    %Calculate interspike intervals (i.e., the log difference between spiketimes) across motifs
%     s = RA_spiketimes;
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
% binSize = str2double(get(handles.edit_psthBin, 'String')); %3ms by default
edges = 0:binSize:motifLength;

%Determine which cells are selected for analysis
% indx = get(handles.list_cellSelect,'Value');
% selectCells = handles.cellList(indx);

%Calculate and plot PSTH per cell
% cellId = getFieldVectorCell(handles.filtSync, 'cellname');
% for i = 1:numel(selectCells)
%     %Generate per-cell logical mask
%     mask = strcmp(selectCells{i}, cellId);
%     
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

%RA_spikes

%Initialize variables
spikeCorr = [];
fs = 44150;

%Set correlation constants
gaussWidth=0.008; %in seconds
% motifLength = size(handles.template,3)/1000; %in seconds
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

