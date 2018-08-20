function [ClusterSet] = ExtractOPTICSClusters(OrderedFiles, SDASet, SUASet, Chi, MinPnts)
%Input:
% OrderedFiles - the ordered-list output of the OPTICS algortihm
%               Structure: [FileIndex CoreDistance ReachabilityDistance]
% SDASet - Set of steep down area start (:,1) and end (:,2) points in OrderedFiles (w.r.t. Chi)
% SUASet - Set of steep up area start (:,1) and end (:,2) points in OrderedFiles (w.r.t. Chi)
% Chi - thresholding value for start and end of steep area (out of 100)
% MinPnts - number of objects in a neighborhood of an object 
%           (minimal number of objects considered as a cluster)
%--------------------
%Output:
% ClusterSet - Set of cluster start (:,1) and end (:,2) points in OrderedFiles (w.r.t. Chi)

%Convert Chi to fraction of 1
Chi = Chi/100;

%Declare output variables/preallocate those of known length
ClusterSet = [];
DownClusterCand = [];
Dmib = [];

%Initialize Index and Mib (Max in between values)
Index = 1; mib = 0;

%Extracts clusters from OrderedFiles based on the presence of corresponding
%Chi-steep down and up regions provided in SDASetand SUASet
while Index<length(OrderedFiles)
    mib = max(mib, OrderedFiles(Index, 3)); %Calculate current point's mib value
    
    %Check Steep Down Candidates
    if ismember(Index, SDASet(:,1))
        Threshold = (mib*(1-Chi)); %Determine threshold for passing SDA point
        SDArow = find(SDASet(:,1)==Index); %Find the row where this index sits
        D = SDASet(SDArow,:); %Pull out the area of interest
        Locmib =mean(OrderedFiles(D(1,1):D(1,2),3));
        
        %Filter out Down Areas that don't meet threshold; include current
        %area in "passing set"
        if ~isempty(DownClusterCand)
            ToFilterOut = find((OrderedFiles(DownClusterCand(:,1),3)*(1-Chi))<mib);
            DownClusterCand = sortrows(setxor(DownClusterCand, DownClusterCand(ToFilterOut,:), 'rows'),1);
            DownClusterCand = sortrows([DownClusterCand; D],1);
        else
            DownClusterCand = [D];
        end
        
        Dmib = [Dmib; Index Locmib];
        
        %Update Index and mib values
        Index = D(1,2)+1;
        mib=OrderedFiles(Index, 3);
        
    %Given previously found Down Candidates, which are valid Up Candidates
    elseif ismember(Index, SUASet(:,1))
        Threshold = (mib*(1-Chi)); %Determine threshold for passing SDA point
        SUArow = find(SUASet(:,1)==Index); %Find the row where this index sits
        U = SUASet(SUArow,:);
        if ~isempty(DownClusterCand)
            ToFilterOut = find((OrderedFiles(DownClusterCand(:,1),3)*(1-Chi))<mib);
            DownClusterCand = sortrows(setxor(DownClusterCand, DownClusterCand(ToFilterOut,:), 'rows'),1);
        end
        
        %Update Index and mib values
        Index = U(1,2)+1;
        mib=OrderedFiles(Index, 3);
        
        %Check all DownClusterCand members for valid pairs with up regions
        if ~isempty(DownClusterCand)
        for p=1:size(DownClusterCand, 1)
            IsMinSize = (U(1,2)-DownClusterCand(p,1)>=MinPnts);
            %StartTh = (OrderedFiles(DownClusterCand(p,1),3)*(1-Chi)) <= mib;%Dmib(Dmib(:,1)==DownClusterCand(p,1),2);
            EndTh = (OrderedFiles(U(1,2),3)*(1-Chi)) >= Dmib(Dmib(:,1)==DownClusterCand(p,1),2);
%            JitterFactor = min(10,length(DownClusterCand(p,1):DownClusterCand(p,2)));
%             if EndTh == 0
%                 z=1;
%                 while z<=JitterFactor && EndTh==0
%                     JitterStart = DownClusterCand(p,1)+z;
%                     Locmib =max(OrderedFiles(JitterStart:DownClusterCand(p,2),3));
%                     EndTh = (OrderedFiles(U(1,2),3)*(1-Chi)) >= Locmib;
%                     z=z+1;
%                 end
%             end
             %%"Pair is Valid" Code%%
             
            if IsMinSize && EndTh %&& StartTh 
                if (OrderedFiles(DownClusterCand(p,1),3)*(1-Chi))>=(OrderedFiles(U(1,2)+1,3))
                    [ReachDiff Offset] =min(abs(OrderedFiles(DownClusterCand(p,1):DownClusterCand(p,2),3)-OrderedFiles(U(1,2),3)));
                    Start = DownClusterCand(p,1)+Offset-1;
                    End = U(1,2);
                elseif OrderedFiles(DownClusterCand(p,1),3)<=(OrderedFiles(U(1,2)+1,3)*(1-Chi))
                    [ReachDiff Offset] =min(abs(OrderedFiles(U(1,1):U(1,2),3)-OrderedFiles(DownClusterCand(p,1),3)));
                    Start = DownClusterCand(p,1);
                    End = U(1,1)+Offset-1;
                else
                    Start = DownClusterCand(p,1);
                    End = U(1,2);
                end
                ClusToAdd = [Start End];
                ClusterSet = [ClusterSet; ClusToAdd];
            end
        end
        end
    else
        Index=Index+1;
    end
end

a=1

end %end function