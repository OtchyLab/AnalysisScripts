function [OrderedFiles, Eps]=OPTICS(DataSet, EpsIn, MinPnts)
% Input: 
% DataSet - data set (m,n); m-objects, n-variables
% EpsIn - user supplied neighborhood radius
%       if not known avoid this parameter or put [] (future plans)
% MinPnts - number of objects in a neighborhood of an object 
%           (minimal number of objects considered as a cluster)
% Output:
% OrderedFile - set of objects ordered by reachability distance
% Structure of OrderedFile is [PermFileIndex, CoreDist, ReachDist]

%Preallocate memory for improved processing speed
Touched = zeros(length(DataSet),1); %file keeps track of which objects have been processed
OrderedFiles = zeros(length(DataSet),3); %file tracks the ordering of objeccts
                                    
%Check Eps value
if isempty(EpsIn)
   Eps = FindEps(DataSet, MinPnts)
else
    Eps = EpsIn;
end

%Sequentially cycles through all objects until all have been "touched"
%Processing (obviously) depends on ordering of objects in DataSet, but the
%ultimate outcome is independent of the DataSet ordering because the entire
%DataSet is sorted for optimal ordering.
a=0;
for PntInd=1:length(DataSet)
    if Touched(PntInd)==0
        a=a+1
        [OrderedFiles, Touched]=ExpandClusterOrder(DataSet, PntInd, Eps, MinPnts, OrderedFiles, Touched);
    end
end


end %end function

function [Eps] = FindEps(x, k)
%Function attempts to estimate an optimal value for Eps from size of
%DataSet and neighborhood members threshold.

[m,n]=size(x);

Eps=((prod(max(x)-min(x))*k*gamma(.5*n+1))/(m*sqrt(pi.^n))).^(1/n);

end %end function

