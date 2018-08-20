function [ClusteredFiles, ClusterSet, Eps]=sortOPTICS(DataSet, EpsIn, MinPnts, Chi)
% Input: 
% DataSet - data set (m,n); m-objects, n-variables
% EpsIn - user supplied neighborhood radius
%       if not known avoid this parameter or put [] (future plans)
% MinPnts - number of objects in a neighborhood of an object 
%           (minimal number of objects considered as a cluster)
% Chi - thresholding value for start and end of steep area (out of 100)
%--------------------------------------------------------------------------
% Output:
% ClusteredFile - set of objects ordered by reachability distance
%       Structure of OrderedFile is [PermFileIndex, CoreDist, ReachDist, ClusterNum]
% ClusterSet - Set of cluster start (:,1) and end (:,2) points in
%       OrderedFiles (w.r.t. Chi)
% Eps - used neighborhood radius
%--------------------------------------------------------------------------

%Run PCA on the imported DataSet to 
[PCA_set] = getPCA(DataSet);

%Run OPTICS sorting algorithm and collect output
[OrderedFiles, Eps]=OPTICS(PCA_set, EpsIn, MinPnts);

%Identify Chi-steep regions from OrderedFiles
[SDASet, SUASet, IndexIsDownUp] = SteepRegions(OrderedFiles, Chi, MinPnts);

%Extract clusters fromt the set of Chi-steep regions
[ClusterSetRaw] = ExtractOPTICSClusters(OrderedFiles, SDASet, SUASet, Chi, MinPnts);

%Filter the found set of clusters by size --> only keep those with smaller
%than 25% of the total number of syllables
ClusDiff = ClusterSetRaw(:,2)-ClusterSetRaw(:,1);
SylFrac = .25;
ClusterSet = ClusterSetRaw(find(ClusDiff<(SylFrac*length(DataSet))),:);

%Create discretized attribute image for the entire dataset (hence, the [])
[AttributeImage] = ExtractAttributePlots(OrderedFiles, DataSet, []);

%Create new figure and plot data
PlotReachAttribClust(OrderedFiles, AttributeImage, ClusterSet)

ClusteredFiles = OrderedFiles;
a=1;

end %function end
