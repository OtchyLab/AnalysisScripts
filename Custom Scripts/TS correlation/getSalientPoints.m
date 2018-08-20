function [pnts,locs] = getSalientPoints(TS)
%Takes in a timeseries (TS) and return a set of points (in x-y space) that gives
%the salient peaks and valleys.  This method assumes that all values of TS
%are positive.

lengthTS = length(TS);
minSpaceDist = 3;
minGradient = 0.005;
checkDist = 3;

%Get Peaks
[peaks,pLocs]=findpeaks(TS,'MINPEAKDISTANCE',minSpaceDist);
pmask = zeros(length(peaks),1);
for i = 1:length(peaks)
    steepLeft = (TS(max(1,pLocs(i)-checkDist)) <= (peaks(i)-minGradient));
    steepRight = (TS(min(pLocs(i)+checkDist,lengthTS))) <= (peaks(i)-minGradient);
    if steepLeft && steepRight
        pmask(i) = 1;
    end
end

%Get Valleys
[valleys,vLocs]=findpeaks(-TS,'MINPEAKDISTANCE',minSpaceDist);
valleys = -valleys;
vmask = zeros(length(valleys),1);
for i = 1:length(valleys)
    steepLeft = (TS(max(1,vLocs(i)-checkDist)) >= (valleys(i)+minGradient));
    steepRight = (TS(min(vLocs(i)+checkDist,lengthTS))) >= (valleys(i)+minGradient);
    if steepLeft && steepRight
        vmask(i) = 1;
    end
end

%Merge peaks and valleys into single set
tempLoc = [pLocs(logical(pmask)),vLocs(logical(vmask))];
tempPnts = [peaks(logical(pmask)),valleys(logical(vmask))];

[locs,indx] = sort(tempLoc);
pnts = tempPnts(indx);

%Plot output...
figure
plot(TS); hold on; axis tight; ylim([0,1])
scatter(pLocs,peaks,'or')
scatter(vLocs,valleys,'og')
scatter(locs,pnts,'+k')



