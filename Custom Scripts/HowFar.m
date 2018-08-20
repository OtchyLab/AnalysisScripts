function [Dshort, Dlong] = HowFar(Sa, Sb)
%S1, S2 - a single array containing all of the acoustic features
%extracted for every syllable requested.  This array is (likely) the product
%of running SlicerDicer and scaling the output with ScaleFeats. The format
%for these arrays are:
%       S*(:,1) - Duration
%       S*(:,2:5) - Fundamental Frequency for time bins 1-4
%       S*(:,6) - Time to Half-Peak Amplitude
%       S*(:,7:10) - Frequency Slope for time bins 1-4
%       S*(:,11:13) - Amplitude Slope between time bins 1-4
%       S*(:,14:17) - Spectral Entropy for time bins 1-4
%       S*(:,18:21) - Temporal Entropy for frequency bins 1-4
%       S*(:,22:25) - Spectro-temporal Entropy for time-frequency quadrants 1-4
%       S*(:,26:29) - Amplitude Modulation for time bins 1-4
%       S*(:,30:33) - Frequency Modulation for time bins 1-4

%Test input array for formatting
[m, n] = size(Sa);
[o, p] = size(Sb);

if m~=1 || n ~=33
    error('Sa not appropriately formatted.  Must be 1x33.')
end

if o~=1 || p ~=33
    error('Sb not appropriately formatted.  Must be 1x33.')
end

%Preallocate variables
%Qa = zeros(4,8);
%Qb = zeros(4,8);
Qa = [];
Qb = [];
Dshort = zeros(4,4);
DlongMat = zeros(4,1);
Dlong = [];

%Format input array into quartile matrices.  For each row of the matrix, format is:
%   Q*(:,1) - Duration
%   Q*(:,2) - Fundamental Frequency
%   Q*(:,3) - Frequency Slope
%   Q*(:,4) - Spectral Entropy
%   Q*(:,5) - Temporal Entropy
%   Q*(:,6) - Spectro-temporal Entropy
%   Q*(:,7) - Amplitude Modulation
%   Q*(:,8) - Frequency Modulation
for i = 1:4;
    Qa(i,:) = [Sa(1), Sa(i+1), Sa(i+6), Sa(i+13), Sa(i+17), Sa(i+21), Sa(i+25), Sa(i+29)];
    Qb(i,:) = [Sb(1), Sb(i+1), Sb(i+6), Sb(i+13), Sb(i+17), Sb(i+21), Sb(i+25), Sb(i+29)];
end
%Populate the short distance matrix base on Euclidean dist between
%quartiles
for i=1:4
    for j=1:4
        Dshort(i,j) = sqrt((1/8)*(sum((Qa(i,:)-Qb(j,:)).^2)));
    end
end

%Populate the long distance matrix base on Euclidean dist between
%whole syllables
for i = 1:4;
    DlongMat(i) = Dshort(i,i).^2;
end

%...and create the final scalar value for Dlong
Dlong = sqrt((1/4)*(sum(DlongMat)));

end