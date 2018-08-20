function [mu_ISI, sig_ISI, mu_numBurst, sig_numBurst] = genHVCparams(age)
% function [mu_ISI, sig_ISI, mu_numBurst, sig_numBurst] = genHVCparams(age)
% Given a specified bird age, this function generates physiologically realistic HVC parameters for 
% simulating a network
%
% created by TMO 2017; last modified 10-09-2017
%
% INPUTS:
%   age:            scalar of the developmental age for which HVC parameters
%                   are returned
% OUTPUTS:
%   mu_ISI:         scalar value of the mu/mean of the ISI distribution (in ms)
%   sig_ISI:        scalar value of the sigma/std of the ISI distribution (in ms)
%   mu_numBurst:    scalar value of the mu/mean of the number of spikes in a burst
%   sig_numBurst:   scalar value of the sigma/std of the number of spikes in a burst


% Debugging values
%age = 65;

% Limit check the age variable
if (age < 60) || (age > 205)
    %That's a problem... something is fucked
    display('Something went wrong with the age input... check and try again.')
    beep
    return
end

% Basic HVC characteristics as pulled from Okubo et al (2015) and Hahnloser et al (2002)
% Key age timepoints (dph): 60-65, 90-130, +200
agePoints = [60, 90, 200];

% Defines distribution of HVC burst ISIs (ms)
HVCisi_m = [mean([4.2707, 3.0033]), 1.63, 1.63]; %estimated from Okubo et al and Hahnloser et al 2002
HVCisi_s = [mean([1.9258, 1.7302])  4.76, 4.76]; %estimated from Okubo et al and Hahnloser et al 2002

% Defines distribution of number of spikes in an HVC burst
HVCnum_m = [mean([4.7, 3.8]), 4.5, 4.5]; %estimated from Okubo et al and Hahnloser et al 2002
HVCnum_s = [mean([2.2, 1.3]), 2, 2]; %estimated from Okubo et al and Hahnloser et al 2002

%Based on the passed in age, calculate the distribution parameters via linear interpolation
m_ISI = interp1(agePoints, HVCisi_m, age);
s_ISI = interp1(agePoints, HVCisi_s, age);
m_numBurst = interp1(agePoints, HVCnum_m, age);
s_numBurst = interp1(agePoints, HVCnum_s, age);

%Generate lognormal mu and sigma from the HVC descriptive statistics
mu_ISI = 2 * log(m_ISI) - 0.5 * log(s_ISI + m_ISI^2);
sig_ISI = sqrt(log(s_ISI + m_ISI^2) - 2 * log(m_ISI));
mu_numBurst = 2 * log(m_numBurst) - 0.5 * log(s_numBurst + m_numBurst^2);
sig_numBurst = sqrt(log(s_numBurst + m_numBurst^2) - 2 * log(m_numBurst));

function [prune_vec, mean_vec, std_vec, ages, inh_vec, HVC_vec, HVCisi_m_vec, HVCisi_s_vec, HVCnum_m_vec, HVCnum_s_vec] = expVects
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

%Sparse datapoints for d60, d110, and d200
sp_ages = [42, 62, 90, 200];

%HVC characteristics
sp_HVC = [6, 4, 2, 2]; %age-dependent decrease in HVC ISI (ms)
%sp_HVC = [2, 2, 2, 2]; %age-dependent decrease in HVC ISI (ms)
% sp_HVCisi_m = 10.^[0.7626, 0.6305, 0.4776, 0.4776]; %estimated from Okubo et al 2015 data (log units)
% sp_HVCisi_s = 10.^[0.2765, 0.2846, 0.2381, 0.2381]; %estimated from Okubo et al 2015 data (log units)
sp_HVCisi_m = [5.7890, mean([4.2707, 3.0033]), 1.63, 1.63]; %estimated from Okubo et al and Hahnloser et al 2002
sp_HVCisi_s = [1.8902, mean([1.9258, 1.7302])  4.76, 4.76]; %estimated from Okubo et al and Hahnloser et al 2002

% sp_HVCnum_m = [3.5, 4.7, 3.8, 3.8]; %estimated from Okubo et al 2015 data
% sp_HVCnum_s = [0.57, 2.2, 1.3, 1.3]; %estimated from Okubo et al 2015 data
sp_HVCnum_m = [3.5, mean([4.7, 3.8]), 4.5, 4.5]; %estimated from Okubo et al and Hahnloser et al 2002
sp_HVCnum_s = [0.6, mean([2.2, 1.3]), 2, 2]; %estimated from Okubo et al and Hahnloser et al 2002

%Inhibition scaling
sp_inh = [1, 1, 1, 1]; % age-dependent inhibition based on Olveczky 2011
%sp_inh = [1, 1, 0.625, 0.625]; % age-dependent inhibition based on Sakaguchi 1996

%Tutored
sp_Mean = [30, 50, 75, 80]./1000;
sp_Std  = [22, 30, 60, 44]./1000;
sp_rawPrun = [18, 26.4, 10.8, 10];

%Isolate
% sp_Mean = [29.9, 45.4, 73, 67]./1000;
% sp_Std  = [20.1, 30.4, 72.8, 57.1]./1000;
% sp_rawPrun = [18, 31, 37, 23];

%Linear scaling function for the pruning
%The basic premise is that the value of rho should be proportional to the extent of synapse pruning such that
%min inputs => rho==1   all synapses pruned
%max inputs => rho==0   no synapses pruned
in_max = 50;
in_min = 0;
sp_Prun = ((in_max-in_min)-sp_rawPrun)./in_max; 

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

r1 = linspace(sp_HVCisi_m(1), sp_HVCisi_m(2), s1); r2 = linspace(sp_HVCisi_m(2), sp_HVCisi_m(3), s2); r3 = linspace(sp_HVCisi_m(3), sp_HVCisi_m(4), s3);
HVCisi_m_vec = [r1, r2(2:end), r3(2:end)];

r1 = linspace(sp_HVCisi_s(1), sp_HVCisi_s(2), s1); r2 = linspace(sp_HVCisi_s(2), sp_HVCisi_s(3), s2); r3 = linspace(sp_HVCisi_s(3), sp_HVCisi_s(4), s3);
HVCisi_s_vec = [r1, r2(2:end), r3(2:end)];

r1 = linspace(sp_HVCnum_m(1), sp_HVCnum_m(2), s1); r2 = linspace(sp_HVCnum_m(2), sp_HVCnum_m(3), s2); r3 = linspace(sp_HVCnum_m(3), sp_HVCnum_m(4), s3);
HVCnum_m_vec = [r1, r2(2:end), r3(2:end)];

r1 = linspace(sp_HVCnum_s(1), sp_HVCnum_s(2), s1); r2 = linspace(sp_HVCnum_s(2), sp_HVCnum_s(3), s2); r3 = linspace(sp_HVCnum_s(3), sp_HVCnum_s(4), s3);
HVCnum_s_vec = [r1, r2(2:end), r3(2:end)];