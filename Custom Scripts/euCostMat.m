function C = euCostMat(m,n)
%Calculates the euclidean distance between evenly spaced points in a two-dimensional grid
%To avoid/minimize memory problems with high dimensional vectors.
%Inputs:
%m = vectors 1
%m = vectors 2

%Generate list of all meshpoints
ms = 1:m; ns = 1:n;
x1 = repmat(ms',[n,1]);
x2 = repmat(ns,[m,1]);
x2 = x2(:);
pnts = [x1, x2];

%Calculate the distance between all points
pd = pdist(pnts);

%And reshape to the square matrix format
C = squareform(pdist(pnts));