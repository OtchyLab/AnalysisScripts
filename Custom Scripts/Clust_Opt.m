function [optimal_groups]=Clust_Opt(DataSet,max_clusters,songs);
%PCA'ed DataSet is Passed from the calling program and the optimized number
%of clusters is returned.  The main idea is that this will do iterative
%clustering until an F-test is satisfied.

%[n p] = size(DataSet);
g = 1:max_clusters; %sets the bounds of testing


IDX = []; IDX1 = [];
C = []; C1 = [];
sumd = []; sumd1 = []; 
centroid_dist = [];
edge_dist = [];
small_clust = [];

dist=squareform(pdist(DataSet)); %Calc distance between clusters in PCA-space

opts = statset('Display', 'off'); %set display aparamters for clustering report

for i=2:length(g);
    IDX = genvarname('IDX', who);
    C = genvarname('C', who);
    sumd = genvarname('sumd', who);
    %eval(['[' IDX ',' C ',' sumd '] = kmeans(DataSet, g(' i '), "distance", "sqEuclidean", "emptyaction", "drop", "onlinephase", "on", "options", opts, "replicates", 5, "start", "cluster");']);
    eval(['[' IDX ',' C ',' sumd '] = kmeans(DataSet, g(i), ''distance'', ''sqEuclidean'', ''emptyaction'', ''drop'', ''onlinephase'', ''on'', ''options'', opts, ''replicates'', 5, ''start'', ''cluster'');']);
    
    %Calculate sum of sqrd Euclidean dist between all centroids
    centroid_matrix = eval(['squareform((pdist(' C ')).^2);']);
    centroid_distb = sum(centroid_matrix)/(2*(length(centroid_matrix)-1));
    
    %Calc sqrd Euclidean dist from the edge of each cluster to the edge of
    %every other cluster and the sum of the syllable to syllable distance
    %within in each cluster
    for j = 1:i;
        cur_clust=eval(['find(' IDX '==j);']);
        cur_clust_pnts = DataSet(cur_clust,:);
        if length(cur_clust)>1
            clust_dist(j) = sum(sum(squareform(pdist(cur_clust_pnts))))/(2*(length(cur_clust)-1));
        else
            clust_dist(j) = 0;
        end
        for k = 1:i
            other = eval(['find(' IDX '==k);']);
            edge_matrix(j,k)=(min(min(dist(cur_clust,other)))).^2;
        end
    end
    
    edge_distb = sum(edge_matrix)/(2*(length(edge_matrix)-1));
    
    %Bin the indexing function and calculate a weighting matrix
    binned = eval(['hist(' IDX ',i);']);
    weights = binned./length(DataSet);
    
    %Calc number of syllables/cluster and threashhold against #songs
    small_clust(i) = (length(find(binned<(0.25*songs))))/i;
    
    %Calc the weighted sum of the mean distances between centroids
    centroid_dist(i) = sum(centroid_distb.*weights);
    
    %Caluclate the weighted sum of the mean intra-cluster distances
    intra_clust_dist(i) = sum(clust_dist.*weights);
    
    %Calculated the weighted sum of the mean distance between cluster edges
    edge_dist(i) = sum(edge_distb.*weights);

end

a=1

optimal_groups = '?';

% F = [];
% Fcrit = [];
% 
% for j=1:length(g)-1
%     for k=2:length(g)
%         lsum = sum(eval(['sumd' int2str(j)]));
%         hsum = sum(eval(['sumd' int2str(k)]));
%         numer = (lsum-hsum)/hsum;
%         denom  = ((n-j)/(n-(k)))*(((k)/j)^(2/p))-1;
%         F(j,k) = numer/denom;
%         Fcrit(j,k) = finv(.99, (p*(k-j)), (p*(n-k)));
%     end
% end

a=1


