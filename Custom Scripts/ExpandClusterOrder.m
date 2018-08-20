function [OrderedFile, Touched]=ExpandClusterOrder(DataSet, PntInd, Eps, MinPnts, OrderedFile, Touched)
% Input: 
% DataSet - data set (m,n); m-objects, n-variables
% PntInd - index of current object being processed
% Eps - neighborhood radius
% MinPnts - number of objects in a neighborhood of an object 
%           (minimal number of objects considered as a cluster)
% OrderedFile - set of objects ordered by reachability distance (up to the
%              current interation)
% Touched - set of objects that have been processed
%-----------------------------------------
% Output:
% OrderedFile - set of objects ordered by reachability distance
% Touched - set of objects that have been processed

% Preallocate/declare local variables
D = [];
Neighbors = [];
ReachDist = [];
CoreDist = [];
NeighDist = [];
OrderSeeds = [];

%Find all of the neighbors for the current point
CurrentPnt = DataSet(PntInd,:); %retrieve current object
D = dist(CurrentPnt, DataSet); %calculate euclidean distance to all other points
Neighbors = find(D<=Eps);%find those objects in the current object's neighborhood

%*****Question: Should Neighbors include the current object itself?  Does 
%*****it matter? Figure 4 suggests it should be included.*****

%Mark current point as "processed"
Touched(PntInd) = 1;
ReachDist(PntInd) = NaN;

%Determine Core Distance (if any)
if length(Neighbors)<MinPnts
    CoreDist(PntInd) = NaN;
else
    NeighDist = [Neighbors' D(Neighbors)'];
    OrderedDist = sortrows(NeighDist, 2);
    CoreDist(PntInd) = OrderedDist(MinPnts, 2);
end

%Write current object parameters to first unused position in OrderedFile
InsertPos = find(OrderedFile(:,1) == 0, 1, 'first');
if InsertPos == 1
    OrderedFile = [PntInd, CoreDist(PntInd), ReachDist(PntInd);OrderedFile(2:end, :)];
else
    OrderedFile = [OrderedFile(1:InsertPos-1, :); PntInd, CoreDist(PntInd), ReachDist(PntInd); OrderedFile(InsertPos+1:end, :)];
end
b=0;
%If the current object is not a core object for this value of Eps, return 
%to the calling function (OPTICS.m); if not, proceed to the sorting function
if isnan(CoreDist(PntInd));
    return
else
    [OrderSeeds] = UpdateSeeds(NeighDist(:,1), CoreDist(PntInd), NeighDist(:,2), Touched, OrderSeeds);
    if ~isempty(OrderSeeds)
    while isempty(find(Touched(OrderSeeds(:,1))==0, 1)) ~= 1
        b=b+1;
        if b==1000
            d=1
        elseif b==2000
            d=1
        elseif b==3000
            d=1
        elseif b==4000
            d=1    
        elseif b==5000
            d=1    
        elseif b==6000
            d=1  
        end
        %Set CurrentPnt to the object in OrderSeeds with smallest
        %Reachability Dist that hasn't yet been touched
        Pointer = find(Touched(OrderSeeds(:,1))==0, 1, 'first');
        CurrentPnt = OrderSeeds(Pointer, 1);
        CurrentObj = DataSet(CurrentPnt,:);
        D = dist(CurrentObj, DataSet);
        Neighbors = find(D<=Eps);%find those objects in the current object's neighborhood
        ReachDist(CurrentPnt) = OrderSeeds(Pointer,2);
        Touched(CurrentPnt) = 1;
        
        %Determine Core Distance of the new point(if any)
        if length(Neighbors)<MinPnts
            CoreDist(CurrentPnt) = NaN;
        else
            NeighDist = [Neighbors' D(Neighbors)'];
            OrderedDist = sortrows(NeighDist, 2);
            CoreDist(CurrentPnt) = OrderedDist(MinPnts, 2);
        end
        
        %Write current object parameters to first unused position in OrderedFile
        InsertPos = find(OrderedFile(:,1) == 0, 1, 'first');
        if InsertPos == 1
            OrderedFile = [PntInd, CoreDist(PntInd), ReachDist(PntInd); OrderedFile(2:end, :)];
        else
            OrderedFile = [OrderedFile(1:InsertPos-1, :); CurrentPnt, CoreDist(CurrentPnt), ReachDist(CurrentPnt); OrderedFile(InsertPos+1:end, :)];
        end
        
        %If the currentPnt is a core object, add these points to the
        %OrderSeeds List and reshuffle it according to reachability
        if isnan(CoreDist(CurrentPnt)) ~= 1
            [OrderSeeds] = UpdateSeeds(NeighDist(:,1), CoreDist(CurrentPnt), NeighDist(:,2), Touched, OrderSeeds);
        end
    end
    end
end
end %end function


function [D]=dist(i,x)
% Ultrafast Euclidean Distance function.
% Calculates distance between some point i in n-dimensional space and every
% point defined in dataset x.
[m,n]=size(x);
D=sqrt(sum((((ones(m,1)*i)-x).^2)'));
if n==1
   D=abs((ones(m,1)*i-x))';
end
end