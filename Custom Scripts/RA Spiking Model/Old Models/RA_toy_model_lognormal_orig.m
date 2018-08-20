clear all;

dt = 0.2;

raster = 1;

t_stop = 1000;  % length of each trial (ms);
N_HVC = 100; % Number of HVC neurons
LMAN_rate = 80; % Rate of LMAN input (Hz)
t = dt : dt : t_stop;

max_mean = 75/1000;
min_mean = 50/1000;

max_std = 80 /1000;
min_std = 30/1000;

min_prune = 0.1;
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

prune_vec = [prune_vec_head , prune_vec, prune_vec_tail];
mean_vec = [mean_vec_head, mean_vec, mean_vec_tail];
std_vec = [std_vec_head, std_vec, std_vec_tail];


tau_s = 5;      % Time constant of HVC-RA and LMAN-RA synapses (ms)
tau_m = 20;     % RA synaptic time constant (ms)
V_th = 20;      % RA firing threshold (mV)
V_r = 0;
t_ref = 1;    % Refractory period of RA (ms)
tau_NMDA = 100; %ms;
c_Mg = 0.7;
LMAN_res = 2* 130;
HVC_res = 2* 130;


noise_std = 0;
W_LMAN =  150/1000;   % LMAN-RA synaptic strength (mV) For NMDA
inh_LMAN = 0;
inh_HVC_ratio = 800;
N_ratio = 0.9; % NMDA/(NMDA+AMP) ratio


% noise_std = 0;
% W_LMAN = 90;   % LMAN-RA synaptic strength (mV) For NMDA
% inh = 150;
% N_ratio = 0.85; % NMDA/(NMDA+AMP) ratio


% noise_std = 0;
% W_LMAN = 220;   % LMAN-RA synaptic strength (mV) for AMPA
% inh = 45;
% N_ratio = 0; % NMDA/(NMDA+AMP) ratio(


% noise_std = 150;
% inh = 100;
% W_LMAN = 0;   % LMAN-RA synaptic strength (mV) for AMPA
% N_ratio = 0; % NMDA/(NMDA+AMP) ratio




pool = 1;

W_HVC_mat = zeros(length(mean_vec),pool*N_HVC);

%pruning_vec = 2/3;

CC_vec =  zeros(size(prune_vec));
CC_vec_2 = zeros(size(prune_vec));
rate_vec = zeros(size(prune_vec));
rec = zeros(size(t));

for pp = 1 : pool
    %clc
    if (raster == 1)
        figure(2);
        clf
        figure(3);
        clf
    end;
    for p = 1 : length(prune_vec)

        

        pruning = prune_vec(p);
        mean_w = mean_vec(p);
        var_w = std_vec(p)^2;
        
        mu = 2 * log(mean_w) - 0.5 * log(var_w + mean_w^2);
        sig = sqrt(log(var_w + mean_w^2) - 2 * log(mean_w));
        
        
        
        trials = 100; % number of trials
        PSTH_bin_size = 20;
        
        
        
        
        % one burst of HVC
        HVC_burst_template = repmat([1, zeros(1, round((2-dt)/dt))], 1 , 5);
        
        
        
        rec_c = zeros(1,length(t) * trials);
        
        % HVC-RA synaptic weights
        %W_HVC = exp(randn(1,N_HVC) * sig + mu);
        W_HVC = lognrnd(mu, sig,1,N_HVC);
        disp([p, mean(W_HVC(W_HVC>0)), std(W_HVC(W_HVC>0))]);
        % Pruning
        W_HVC = W_HVC .* (rand(size(W_HVC)) >= pruning);
        
        inh_HVC = mean(W_HVC) * inh_HVC_ratio;
        inh = inh_LMAN + inh_HVC;
        
        
        % Geberating the whole HVC spikes (modulated by the corresponding rate)
        HVC_spikes = [];
        for i = 1 : N_HVC
            HVC_spikes = [ HVC_spikes, HVC_burst_template * W_HVC(i)];
        end;
        % Equalizing the length to t
        HVC_spikes = [HVC_spikes, zeros(size(t))];
        HVC_spikes = HVC_spikes(1: length(t));
        RA_spikes = zeros(trials, length(t));
        
        %%
%         HVC_spikes = zeros(size(HVC_spikes));
%         HVC_spikes(10) = 0;
%         W_HVC = ones(size(W_HVC));
%         W_LMAN = 1;   % LMAN-RA synaptic strength (mV) For NMDA
%         inh = 0;
%         N_ratio = 1;
        %%
        prev_spike_time = 0;
        ISI = zeros(trials, length(t));
        % loop on trials
        for i = 1 : trials
            
            % LMAN input to RA neuron
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
            end;
            ref = 0;
            
            %% loop on time
            prev_tt = 1;
            for tt = 2 : length(t)
                
                x_Mg = 1/(1 + exp(- 0.062 * (V_RA(tt-1)-70)) * c_Mg/3.57);
                dLMAN_AMPA_input_dt = (- LMAN_AMPA_input(tt - 1) + tau_s * (1-N_ratio) * W_LMAN * LMAN_spikes(tt)/dt)/tau_s;
                dLMAN_NMDA_input_dt = (- LMAN_NMDA_input(tt - 1) + tau_NMDA * N_ratio * x_Mg * W_LMAN * LMAN_spikes(tt)/dt)/tau_NMDA;
                dHVC_input_dt = ( - HVC_input(tt - 1) + tau_s * HVC_spikes(tt)/dt)/tau_s;
                
                LMAN_AMPA_input(tt) = LMAN_AMPA_input(tt - 1) + dt * dLMAN_AMPA_input_dt;
                LMAN_NMDA_input(tt) = LMAN_NMDA_input(tt - 1) + dt * dLMAN_NMDA_input_dt;
                HVC_input(tt) = HVC_input(tt - 1) + dt * dHVC_input_dt;
                
                I_noise = noise_std/sqrt(dt) * randn;
                
                %% RA_neuron
                
                dV_RA_dt = ( - V_RA(tt-1) + LMAN_AMPA_input(tt) * LMAN_res + LMAN_NMDA_input(tt) * LMAN_res + HVC_input(tt) * HVC_res ...
                             - inh + I_noise)/tau_m;
                
                if ref<=0
                    V_RA(tt) = V_RA(tt-1) + dt * dV_RA_dt;
                else
                    V_RA(tt) = 0;
                end;
                ref = ref - dt;
                %% firing condition;
                if (V_RA(tt) >=V_th) && (ref<=0)
                    RA_spikes(i,tt) = 1;
                    V_RA(tt) = V_r;
                    ref = t_ref;
                    ISI(i, prev_tt : tt) = t(tt) - prev_spike_time;
                    prev_spike_time = t(tt);
                    prev_tt = tt;
                end;
                
                
                
                
            rec(tt) = LMAN_NMDA_input(tt);     
            end;
            
            % RA raster plot
            if (i<=100) && (raster == 1) && ((p==1) || (p == (length(prune_vec))))
                RA_spike_times = RA_spikes(i,:) .* t;
                RA_spike_times(RA_spike_times == 0) = [];
                if p == 1
                    figure(2);
                    xlabel(num2str(pruning));
                else
                    figure(3);
                    xlabel(num2str(pruning));
                end;
                hold on;
                X = [RA_spike_times; RA_spike_times];
                Y = [(i-1) * 3 * ones(size(RA_spike_times)); ((i-1) * 3 + 2) * ones(size(RA_spike_times))];
                line(X,Y,'color','k');
                hold on;
%                 plot(RA_spike_times, i * ones(size(RA_spike_times)) ,'marker','.', 'markersize', 7, 'linestyle','none','color','k');
%                  hold on;
%                 %disp('raster')
            end;
            
        end;
        
        bin_num = round(PSTH_bin_size/dt);
        
        for i = 1 : trials
            RA_spikes_bined(i,:) = sum(reshape(RA_spikes(i,:), bin_num, round(length(t)/bin_num)),1);
        end;
        
        
        
        
        
        
        
        
        RA_rate_mean = 1000* mean(RA_spikes_bined,1)/PSTH_bin_size;
        RA_rate_std =  std(1000* RA_spikes_bined /(PSTH_bin_size),1)/sqrt(PSTH_bin_size);
        
        % RA_rate_mean = 1000 * reshape(sum(RA_spikes,1), bin_num, round(length(t)/bin_num))/ (trials * dt);
        % RA_rate_mean = mean(RA_rate_mean,1);
        
        
        RA_rate_smooth = smoothts(1000./ISI,'g',1000,10/dt);
        RA_rate_smooth = 1000./ISI;
        RA_rate_smooth(isnan(RA_rate_smooth)) = 0;
        RA_rate_smooth(isinf(RA_rate_smooth)) = 0;
        
        C = corrcoef(RA_rate_smooth');
        
        
        %C = corrcoef(RA_spikes_bined');
        C = C .* (1 - eye(trials));
        C(isnan(C)) = 0;
        CC = sum(C(:)) / (trials^2 - trials);
%         if RA_rate_mean <20
%             CC = CC_vec(p);
%         end;
        CC_vec(p) =  ((pp - 1) * CC_vec(p) + CC)/pp;
        CC_vec_2(p)= ((pp - 1) * CC_vec_2(p) + CC^2)/pp;
        rate_vec(p) = ((pp - 1) * rate_vec(p) + mean(RA_rate_mean))/pp;
        %[mean(LMAN_AMPA_input), std(LMAN_AMPA_input), 1000*sum(sum(RA_spikes))/t_stop/trials]
%         W_HVC_mat(p,( (pp-1) * N_HVC + 1): pp * N_HVC) = W_HVC;
%         W = W_HVC_mat(p,:);
%         figure(4);
%         hand = subplot(length(prune_vec),1,p);
%         hist(W(W>0),20);
%         set(hand,'xlim', [0, 0.3]);
%         disp([p, mean(W(W>0)), std(W(W>0))]);          
    end;
    

figure(1);

plot(prune_vec, rate_vec,'r');
hold on;
plot(prune_vec, CC_vec*100,'b');
hold off;
drawnow

end;

