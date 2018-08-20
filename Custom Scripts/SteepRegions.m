function [SDASet SUASet IndexIsDownUp] = SteepRegions(OrderedFiles, Chi, MinPnts)
%Input:
% OrderedFiles - the ordered-list output of the OPTICS algortihm
%               Structure: [FileIndex CoreDistance ReachabilityDistance]
% Chi - thresholding value for start and end of steep area (out of 100)
% MinPnts - number of objects in a neighborhood of an object 
%           (minimal number of objects considered as a cluster)
%-----------------
%Output:
% SDASet - Set of steep down area start (:,1) and end (:,2) points in OrderedFiles (w.r.t. Chi)
% SUASet - Set of steep up area start (:,1) and end (:,2) points in OrderedFiles (w.r.t. Chi)
% IndexIsDownUp - Index of up/down points for troubleshooting

%Convert Chi to fraction of 1
Chi = Chi/100;

%Declare output variables/preallocate those of known length
SDASet = []; 
SUASet = [];
IndexIsDownUp = zeros(length(OrderedFiles),2);

%Calculate Up/Down Pos for each point in OrderedFiles
%Stucture is logical = [IsaDownPoint IsaUpPoint]
i=1:length(OrderedFiles)-1;
IndexIsDownUp(i,1) = (OrderedFiles(i,3)*(1-Chi)) >= OrderedFiles(i+1,3); %Point is Chi-steep down
IndexIsDownUp(i,2) = OrderedFiles(i,3) <= (OrderedFiles(i+1,3)*(1-Chi)); %Point is Chi-steep up
IndexIsDownUp(i,3) = OrderedFiles(i,3);

%Determine the breakup into Steep UP Regions
SetUpPoints = find(IndexIsDownUp(:,2)==1);
j=1;
while j<=length(SetUpPoints)-1
    ScanRange = SetUpPoints(j)+1:SetUpPoints(j+1); %Candidate Chi-steep Upward range selected
    ind = 1:length(ScanRange);
    AtLeastAsHigh = zeros(length(ScanRange),1);
    AtLeastAsHigh(ind)= OrderedFiles(ScanRange,3)>=OrderedFiles(ScanRange-1,3); %Check to assure that thhis range contains no down turns
    RangeContPass = isequal(sum(AtLeastAsHigh),length(AtLeastAsHigh)); %if all are upward facing, set RangeContPass
    
    if SetUpPoints(j+1)-SetUpPoints(j)>MinPnts
        NotUpsRange = SetUpPoints(j):SetUpPoints(j+1)-MinPnts; %Candidate Chi-steep Upward range selected
        NotUps = zeros(length(NotUpsRange),1);
        for ind = 1:length(NotUpsRange);
        NotUps(ind) = sum(IndexIsDownUp(NotUpsRange(ind):NotUpsRange(ind)+MinPnts-1,2))==0;
        end
        RangeNotUpsPass = sum(NotUps)==0; %if there are no MinPnts-length static regions, set RangeNotUpsPass
    else
        RangeNotUpsPass = 1; %if there are no MinPnts-length static regions, set RangeNotUpsPass
    end
    
    %If all of these conditions are met, the set quaifies as a Chi-Steep Up Area
    if RangeContPass==1 && RangeNotUpsPass==1
         if isempty(SUASet)
            SUASet = [SetUpPoints(j) SetUpPoints(j+1)];
         else
            [m, n] = size(SUASet);
            if isequal(SUASet(m,n), SetUpPoints(j))
                SUASet(m,n) = SetUpPoints(j+1);
            else
                SUASet = [SUASet; SetUpPoints(j) SetUpPoints(j+1)];
            end
         end
    end 
    j=j+1;
end

%Determine the breakup into Steep DOWN Regions
SetDownPoints = find(IndexIsDownUp(:,1)==1);
j=1;
while j<=length(SetDownPoints)-1
    ScanRange = SetDownPoints(j)+1:SetDownPoints(j+1); %Candidate Chi-steep Downward range selected
    ind = 1:length(ScanRange);
    AtLeastAsHigh = zeros(length(ScanRange),1);
    AtLeastAsHigh(ind)= OrderedFiles(ScanRange,3)<=OrderedFiles(ScanRange-1,3); %Check to assure that thhis range contains no up turns
    RangeContPass = isequal(sum(AtLeastAsHigh),length(AtLeastAsHigh)); %if all are downward facing, set RangeContPass
    
    if SetDownPoints(j+1)-SetDownPoints(j)>MinPnts
        NotDownsRange = SetDownPoints(j):SetDownPoints(j+1)-MinPnts; %Candidate Chi-steep Downward range selected
        NotDowns = zeros(length(NotDownsRange),1);
        for ind = 1:length(NotDownsRange);
        NotDowns(ind) = sum(IndexIsDownUp(NotDownsRange(ind):NotDownsRange(ind)+MinPnts-1,1))==0;
        end
        RangeNotDownsPass = sum(NotDowns)==0; %if there are no MinPnts-length static regions, set RangeNotDownsPass
    else
        NotDownsRange = SetDownPoints(j):SetDownPoints(j+1); %Candidate Chi-steep Upward range selected
        RangeNotDownsPass = 1; %if there are no MinPnts-length static regions, set RangeNotDownsPass
    end
    
    %If all of these conditions are met, the set quaifies as a Chi-Steep Up Area
    if RangeContPass==1 && RangeNotDownsPass==1
         if isempty(SDASet)
            SDASet = [SetDownPoints(j) SetDownPoints(j+1)];
         else
            [m, n] = size(SDASet);
            if isequal(SDASet(m,n), SetDownPoints(j))
                SDASet(m,n) = SetDownPoints(j+1);
            else
                SDASet = [SDASet; SetDownPoints(j) SetDownPoints(j+1)];
            end
         end
    end 
    j=j+1;
end

% for k=1:length(SDASet)
% RangetoSet=SDASet(k,1):SDASet(k,2);
% IndexIsDownUp(RangetoSet,4)=666;
% end
% for k=1:length(SUASet)
% RangetoSet=SUASet(k,1):SUASet(k,2);
% IndexIsDownUp(RangetoSet,5)=777;
% end
  
end   % temp end function
