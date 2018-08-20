function SampleClusters(OrderedFiles, SylParam, ClusterSet, n)

%Ask user which file to find the dated folders of individuated syllables
parentsource = uigetdir(pwd,'Source Folder');
%parentsource = 'F:\Rd278 Stripped Syllables';

%Ask user where to store selected syllables
parentdest = uigetdir(pwd,'Destination Folder');

if n>100
    n=100;
end

SylNames = ['A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'; 'K'; 'L'; 'M'; 'O'; 'P'];
num_clus = length(ClusterSet);
IndicesUsed=zeros(num_clus,n);

for i=1:num_clus
    SampleRange = ClusterSet(i,1):ClusterSet(i,2);
    f = ceil((ClusterSet(i,2)-ClusterSet(i,1)).*rand(100,1))+ClusterSet(i,1);
    IndicesUsed(i,:)=f(1:n)';
    SamplePull = OrderedFiles(f(1:n),1);
    for j=1:n
        source = [parentsource '\' eval(['SylParam(' num2str(SamplePull(j)) ').name'])];
        filenew = eval(['SylParam(' num2str(SamplePull(j)) ').name']);
        filenew = regexprep(filenew, '.wav', ['_' SylNames(i) '.wav']);
        dest_app = [parentdest '\' SylNames(i)];
        if ~isdir(dest_app)
            mkdir(dest_app);
        end
        dest_app = [dest_app '\' filenew];
        copyfile(source, dest_app);
%         eval(['FeatureSet_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '_' SylNames(i) '(j,:) = DataSet_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '(:,' num2str(index(j)) ');']);
%         eval(['SylSet_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '_' SylNames(k) '(j,:) =  filenew;']);
        cd(parentdest);
%         save(['SummaryData_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '_' SylNames(k) '.mat'], ['FeatureSet_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '_' SylNames(k)], ['SylSet_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '_' SylNames(k)]);  
    end
end

end %end function