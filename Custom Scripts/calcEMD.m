function emd = calcEMD(Pprob,Qprob,D)
%Calculate the Earth Mover Distance between two probability distributions (P and Q) using Ofir Pele's C-code
%found at http://www.ariel.ac.il/sites/ofirpele/FASTEMD/code/ . For now, this is exclusive written for the 
%heatmaps for the inactivations project

if nargin < 3 %i.e., Distance matrix not passed
%     load('CostMat-thresh40_halfBins.mat'); %load predefined distance matrix for maxTime: 300ms, 2.5ms bins
%     load('CostMatrix 281x161.mat'); %load predefined distance matrix for maxTime: 400ms, 2.5ms bins
end

%Number of fake samples
numSamp = 10000;

%Convert probability distributions to same-size histograms
P = round(Pprob*numSamp);
Q = round(Qprob*numSamp);

%Check the sizing of the histos to assure they are the same
R= size(P,1);
C= size(P,2);
if (~(size(Q,1)==R&&size(Q,2)==C))
    error('Size of images should be the same');
end

%%This is the code for generating the distance/cost matrix -- not used speed/processing reasons, but here for reference
% maxFact = floor(min([R, C])/2);
% COST_MULT_FACTOR= 1;
% THRESHOLD= maxFact*COST_MULT_FACTOR;
% D= zeros(R*C,R*C,'int32');
% j= 0;
% for c1=1:C
%     for r1=1:R
%         j= j+1;
%         i= 0;
%         for c2=1:C
%             for r2=1:R
%                 i= i+1;
%                 D(i,j)= min( [THRESHOLD (COST_MULT_FACTOR*sqrt((r1-r2)^2+(c1-c2)^2))] );
%             end
%         end
%     end
% end

%Mass penalty (-1 == max D penalty)
extra_mass_penalty= int32(-1);

% flowType= int32(3);

%Convert the histograms to the vectors (consistent with the distance matrix) 
P= int32(P(:));
Q= int32(Q(:));

%Calculate the EMD using the C-mex scripts
emd = emd_hat_gd_metric_mex(P,Q,D,extra_mass_penalty);

%Convert to double (from int32) for ease of future use and normalize by the number of samples
emd = double(emd)/mean([sum(P(:)), sum(Q(:))]);



