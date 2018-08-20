function [smOut, res] = sgts(histo)
%This function implements Good-Turing Smoothing as described by Sampson and Gale (1995) "Good-Turing
% Smoothing without Tears" and is based on the psuedocode provided in Appendix II. The follows closely
% Mike Hammond's code at http://dingo.sbs.arizona.edu/~hammond/ling439-f13/gtanal.m
%
%The script takes as input a matrix (or row vector) containing counts at each of the positions.
%Note that  this smoothing is agnostic as to what features/categories the vector/row positions
%correspond to -- it simply cares about the number of tokens for each possible type and all types
%are treated as unique.
%
% [smOut, res] = sgts(histo)
%       histo = the matrix or row vector containing the token counts for each possible type.
%       smOut = the smoothed version of the input histo
%       res = lookup table for converting the original rank values for smoothed estimates
% Timothy Otchy, 4/1/2015

%Takes in input data and converts to rank and frequency of rank vectors (excluding the zero counts!!!)
[xnr,xr] = hist2counts(histo);

%Calculate from the new vectors the total number of tokens in the data set.
xN = sum(xr .* xnr);

% #make averaging transform
%Estimate the frequencies of rare types (i.e., when nr is small) by averaging the surrounding zero 
%and non-zero frequencies. This is the Z-transform described in the original paper.
xnrz = nrzest(xr,xnr);

% Calculate the Linear Good-Turing estimate by linear regression on the log-log plot of Z-transformed seen-type frquencies. 
% (OLD: coefficients in opposite order in matlab)
xcoef = polyfit(log(xr),log(xnrz),1);

% Apply the LGT estimate to the rank vector to generate r*...
xrst = rstest(xr,xcoef);
xrstrel = xrst ./ xr;

% Calculate the Turing estimate
xrtry = (xr == vertcat(xr(2:end)-1,0));

xrstarel = zeros(length(xr),1);

xnrTmp = vertcat(xnr(2:end),0);
xrstarel(xrtry) = (xr(xrtry)+1)./xr(xrtry) .* ...
    xnrTmp(xrtry) ./xnr(xrtry);

% #make switch from Turing to LGT estimates
tursd = ones(length(xr),1);
for i = 1:length(xr)
    if xrtry(i)
        tursd(i) = (i+1)/...
            xnr(i)*sqrt(xnr(i+1)*(1+xnr(i+1)/xnr(i)));
    end
end

%Step through each position in the rank vector and select/apply the appropriate smoothing
xrstcmbrel = zeros(length(xr),1); %Pre-allocate the output vector
useturing = true; %Set initial flag
for r = 1:length(xr)
    %This is the code for switching between using the simple Turing estimate (appropriate for low-frequency counts) and the
    %Good-Turing estimate for higher frequencies. Note that the algorithm can only switch estimators once. 
    if ~useturing
        xrstcmbrel(r) = xrstrel(r);
    elseif abs(xrstrel(r) - xrstarel(r)) * r/tursd(r) > 1.65
        xrstcmbrel(r) = xrstarel(r);
    else
        useturing = false;
        xrstcmbrel(r) = xrstrel(r);
    end
end

% Renormalize the probabilities for observed objects
sumpraw = sum(xrstcmbrel .* xr .* xnr ./ xN);
xrstcmbrel = xrstcmbrel .* (1 - xnr(1)./xN)./sumpraw;

%Generate the lookup table for the smoothed counts.
res = zeros(length(xr)+2,2);
res(1,1) = xN;
res(1,2) = sum(xnr);
res(2,2) = xnr(1)/xN; %This is the reserved frequency count for unseen tokens
res(3:end,1) = xr;
res(3:end,2) = xr .* xrstcmbrel;

%Smooth the original input using the res lookup table
smOut = applySGTS(histo, res);

%Z-transform script
function [res] = nrzest(r,nr)
%Estimate the frequencies for unseen types by taking the linear average of the nearest non-zero neighbors.
%Note that at small r, differences from non-zero nr are likely to be near unity; at large r, differences can be huge.
d = vertcat(1,-1*(r(1:end-1) - r(2:end)));
dr = vertcat((.5 * (d(2:end) + d(1:end-1))),d(end));
res = nr./dr;

%Remap the rank values using the linear regression
function [res] = rstest(r,coef)
res = r .* (1 + 1./r).^(1 + coef(1));

%Convert histogram to rank and counts vectors
function [xnr,xr] = hist2counts(histo)
%Rank vector is simply the integers from 1 to max histogram counts
xr = 1:max(histo(:));

%Generate frequency vector
for i = xr
    xnr(i) = length(find(histo(:) == xr(i)));
end

%Remove zero-counts and transpose
ind = xnr~=0;
xr = xr(ind)';
xnr = xnr(ind)';

function smOut = applySGTS(data, table)
%determine the size of the original dataset
[m, n] = size(data);

%Preallocate output
smOut = zeros(m,n);

%Loop through the entire data and substitute values
for i = 1:size(table,1)
    smOut(data==table(i,1)) = table(i,2); 
end















