function rho = LoCo(TS1,TS2)

if length(TS1)~=length(TS2)
    print('The two time series must have the same length')
    return
end

% Windowing parameters
fs = 44150; % Sampling frequency in Hz
winsizeT = 2; % Window size in ms
winstepT = 1; % Window step in ms

winsize = floor(winsizeT*(fs/1000)); % Convert to sample counts
winstep = floor(winstepT*(fs/1000)); % Convert to sample counts

% Window the two time series
Mat1 = running_windows(TS1,winsize,winstep);
Mat2 = running_windows(TS2,winsize,winstep);

winnum = size(Mat1,1);

% Define variables
gamma1SW = [];
gamma1Exp = [];
gamma2SW = [];
gamma2Exp = [];
beta = .25; %Exponential decay factor; some number between 0 and 1
k = 4; %Number of eigenvectors to use in decomposition

% Calculate autocorrelation matrix for the first window of each series
gamma1SW(1,:,:) = Mat1(1,:)'*Mat1(1,:);
gamma1Exp(1,:,:) = gamma1SW(1,:,:);

gamma2SW(1,:,:) = Mat2(1,:)'*Mat2(1,:);
gamma2Exp(1,:,:) = gamma2SW(1,:,:);

% Calculate local autocovariance matrix for each future window. These are 
% different from the first window b/c I want this value to reflect the
% "neighborhood" that trails the current timepoint.
for i = 2:winnum
    %gamma1SW(1,:,:) = ;
    gamma1Exp(i,:,:) = (beta*squeeze(gamma1Exp(i-1,:,:)))+(Mat1(i,:)'*Mat1(i,:));
    
    %gamma2SW(1,:,:) = ;
    gamma2Exp(i,:,:) = (beta*squeeze(gamma2Exp(i-1,:,:)))+(Mat2(i,:)'*Mat2(i,:));
end

% Calculate the LoCo Score for each of the AutoCov matrices derived above
for i = 1:winnum
    [V,D] = eigs(squeeze(gamma1Exp(i,:,:)),k);
    
end


    