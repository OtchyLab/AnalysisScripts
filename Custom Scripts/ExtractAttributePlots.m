function [AttributeImage] = ExtractAttributePlots(OrderedFiles, DataSet, ClusterSet)
%This function creates a grayscale image that displays the discretized
%attributes for an n-dimensional dataset and aligns that display with the
%reachability plot.
% Input: 
% OrderedFile - set of objects ordered by reachability distance
%        Structure of OrderedFile is [PermFileIndex, CoreDist, ReachDist]
% DataSet - data set (m,n); m-objects, n-variables
% SylParam - structured dataset containing all of the syllable features
%       (both feature-space and z-score space)
% ClusterSet - 2-d array representing the beginning and end index for each
%       identified cluster.  Leave as empty set to calc the entire 
%       OrderedFiles list
% -------------------
% Output:
% AttributeImage - a matrix of 100*n rows and length(OrderedFiles)colums 
%       that represents an image of discretized grayscale values


numAttributes = size(DataSet,2); %Number of different syllable attributes calculated
BarHeight = 2; %Height in pixels for each attribute bar, min=2
SpreadSet = max(DataSet)-min(DataSet); %Calculate the spread for eachof the attributes (here, the z-scores)
MinSet = min(DataSet); %Calculate the minimum value for each attribute (here, the z-scores)
ColorDepth = 255; % Image will be scaled to 8-bit grayscale.  Change to (2^16)-1 this number for 16-bit.

%If ClusterSet is empty, set it as a single cluster
if isempty(ClusterSet)
   ClusterSet = [1, length(OrderedFiles)];
end

%Preallocate the mem space for the image
AttributeImage = zeros(BarHeight*numAttributes, length(OrderedFiles)); 
j = 1:numAttributes; %To vectorize scaling calculation
StartPos(j) = BarHeight*(j-1)+1; %Vertical start position for the AttribImage
EndPos(j) = BarHeight*(j); %Vertical end position for the AttribImage

%Extract the feature attributes from the DataSet
for n=1:size(ClusterSet,1)
    for i = ClusterSet(n,1):ClusterSet(n,2);
        ScaledAttrib(j) = round((ColorDepth)*(DataSet(OrderedFiles(i,1),j)-MinSet(j))./SpreadSet(j)); %Scale features by ColorDepth
        for s = 1:numAttributes
            AttributeImage(StartPos(s):EndPos(s),i) = ScaledAttrib(s); %Insert scaled values along the vertical dimension
        end
    end
end


%imagesc(AttributeImage,[0 ColorDepth]); colormap(gray); %Plot the constructed image using grayscale values


end %end function