function masterHVC
% This code is for creating an age-relative HVC sequence based on the spiking
% characteristics previously described by Okubo et al (2015) and Hahnloser
% (2002). With the debugging parameters uncommented, it can run on it's own
% to generate an HVC spiketimes matrix. (To plot the output, you may need
% to change the boolean flag in the genHVCpool function.
%
% created by TMO 2017; last modified 10-09-2017

% Debugging/test parameters
age = 90; %simulated bird age
nHVC = 100; %number of HVC neurons to simulate
dt = 0.2; %time step (ms)
nodeSpace = 10; %spacing of HVC nodes in the model

%Retrieve the parameters for generating an HVC pool at this age
[mu_ISI, sig_ISI, mu_numBurst, sig_numBurst] = genHVCparams(age);

%Generate the HVC pool from the array
spikeMat = genHVCpool(nHVC, mu_ISI, sig_ISI, mu_numBurst, sig_numBurst, dt, nodeSpace);





