function a=crosscorr1(t1,t2,N_bins);
 
M1 = length(t1);
M2 = length(t2);
D = ones(M2,1)*t1 - t2'*ones(1,M1);
D = D(:);
a=hist(D,N_bins);