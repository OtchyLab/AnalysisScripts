function [OrderSeeds]=UpdateSeeds(Neighbors, CoreDist, NeighborDists, Touched, OrderSeeds)
%The purpose of this function is to create and update a list (OrderSeeds)
%such that it is constantly organized by reachability distance from the
%core object.
% Inputs:
% Neighbors - the DataSet indexes of the Neighbors to the Core Object
% CoreDist - the determined core distance for the core object
% NeighborDist - the Euclidean distance (in n-dimensional space) between
%           `    the core object and each of the neighbors
% Touched - the current list of which objects have been processed
% OrderSeeds - the current list of OrderSeeds
%
% Outputs:
% OrderSeeds - The Reachability-ordered list of Neighbors and distances

%Allocate variable and memo0ry for the lists
ReachDist = NaN(length(Neighbors),1)';


%Cycle through each of the Neighbors to calc core distances and then order
%the list
for i = 1:length(Neighbors)
    if Touched(Neighbors(i))==0 %check to see if the point has been processed yet, if not continue on
        %Candidate Reachability dist is the max of the core object dist or Euclidean distance
        NewReach = max(CoreDist, NeighborDists(i));
        if isempty(OrderSeeds)
            IsNotaMember = 1;
        else
            IsNotaMember=~ismember(Neighbors(i), OrderSeeds(:,1));
        end
        if isnan(ReachDist(i)) && IsNotaMember==1%if you haven't yet set any reachability distance...
            ReachDist(i) = NewReach;% ...this will be the reachability distance
            if isempty(OrderSeeds)
                OrderSeeds = [Neighbors(i), ReachDist(i)];
            else
                InsertPos = find(OrderSeeds(:,2)<=NewReach, 1, 'last');
                if isempty(InsertPos) ||InsertPos == size(OrderSeeds,1)
                    OrderSeeds = [OrderSeeds; Neighbors(i), ReachDist(i)];
                else
                    OrderSeeds = [OrderSeeds(1:InsertPos,:); Neighbors(i), ReachDist(i); OrderSeeds(InsertPos+1:end,:)];
                end
            end
        elseif ismember(Neighbors(i), OrderSeeds(:,1))==1 %If the object is already in the Ordered Seed list...
            %find out where it is
            PlaceInSeeds = find(OrderSeeds(:,1)==Neighbors(i));
            if NewReach < OrderSeeds(PlaceInSeeds,2) %if the current reachDist is less than the stored dist
                %make this the new reachDist
                ReachDist(i) = NewReach; 
                OrderSeeds(PlaceInSeeds,2) = NewReach;
                %Next few lines handle the resorting of OrderSeeds based on the new dist
                ObjToMove = OrderSeeds(PlaceInSeeds,:); 
                TempOrderSeeds = sortrows(setxor(OrderSeeds, ObjToMove, 'rows'),2);
                InsertPos = find(TempOrderSeeds(:,2)<=NewReach, 1, 'last');
                OrderSeeds = [TempOrderSeeds(1:InsertPos,:); ObjToMove; TempOrderSeeds(InsertPos+1:end,:)];
            end
        end
    end
end


end %end function
