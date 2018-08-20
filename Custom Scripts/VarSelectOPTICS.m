function VarSelectOPTICS(end_date, earliest_date, interval)

%Ask user which file to find the dated folders of individuated syllables
%parentsource = uigetdir(pwd,'Source Folder');
parentsource = 'F:\Rd278 Stripped Syllables';

%Ask user where to store selected syllables
%parentdest = uigetdir(pwd,'Destination Folder');
parentdest = 'C:\Documents and Settings\Bird Brain\Desktop\Rd278 Sorted OPTICS';

if isdir(parentdest)
else
    mkdir(parentdest);
end

date = end_date;
SylNames = ['A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'; 'K'; 'L'; 'M'; 'O'; 'P'];

songsource = [parentsource '\' date];

while exist([songsource '\' date '.mat'])~=2
    date = datestr(addtodate(datenum(date), -1, 'day'), 'yyyy-mm-dd');
    songsource = [parentsource '\' date];
end
    
load ([songsource '\' date '.mat']);

while datenum(date, 'yyyy-mm-dd') >= datenum(earliest_date, 'yyyy-mm-dd')
    songsource = [parentsource '\' date];
    if isdir(songsource)
        load ([songsource '\' date '.mat']);
        dest = [parentdest '\' date];

        cur_DataSet = eval(['DataSet_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd')]);
        [coeff score latent] = princomp(cur_DataSet');

        for t=1:length(latent)
            vartot(t)=sum(latent(1:t))/length(latent);
        end
        num_pc = find(vartot>.95, 1); %Determine how many prin components are required for 95% variation
        PCA_set = score(:,1:num_pc); %Reduce working data to only those components
        
        if isdir(dest)
        else
            mkdir(dest);
        end
%         %***************
%         k=8;        %Number of neighbors required to define a cluster
%         Eps = .4;   %Neighborhood radius (leave empty for algorithmic estimate)
%         %***************
%         [cluster, type] = dbscan(PCA_set,k,Eps);
%         clus_set = unique(cluster);
%         num_clus = length(clus_set);

        Eps = 2;
        MinPnts = 150;
        Chi = 0.05;

        [OrderedFiles, Touched]=OPTICS(PCA_set, Eps, MinPnts);
        [SDASet SUASet IndexIsDownUp] = SteepRegions(OrderedFiles, Chi, MinPnts);
        [ClusterSet] = ExtractOPTICSClusters(OrderedFiles, SDASet, SUASet, Chi, MinPnts);
        
        clus_set = 1:length(ClusterSet);
        num_clus = length(clus_set);
        
        cluster = zeros(length(OrderedFiles),1);
        for q=1:num_clus;
            InsertRange = ClusterSet(q,1):ClusterSet(q,2);
            cluster(OrderedFiles(InsertRange,1)) = q;
        end
        
        c = ['r';'b'; 'k'; 'g'; 'y'; 'm'; 'c'];
        mar = ['o'; '+'; 's'];
    

        for u=1:num_clus
            for v=2:num_pc
                figure(v);
                if u<=7
                scatter(PCA_set((cluster==clus_set(u)),1),PCA_set((cluster==clus_set(u)),v),2, c(u), mar(1));
                elseif u>7 & u<=14
                scatter(PCA_set((cluster==clus_set(u)),1),PCA_set((cluster==clus_set(u)),v),3, c(u-7), mar(2));
                else
                scatter(PCA_set((cluster==clus_set(u)),1),PCA_set((cluster==clus_set(u)),v),3, c(u-14), mar(3));
                end
                hold on
                xlabel('PC1');
                ylabel(['PC' num2str(v)]);
                title(['Scatterplot of PC1 vs PC' num2str(v) ' with k=' num2str(k) ' and Eps=' num2str(Eps)]);
            end
        end
        hold off
        
        for k=1:num_clus
            %index = find(0 < simfeat(k,length(test_syls)+1:end) & simfeat(k,length(test_syls)+1:end) <= dist_th);%-length(test_syls);
            index = find(clus_set(k)==cluster); %find all syllables both assigned to the cluster and not an outlier
            sylnum = length(index);
            z_Dur = 0; z_FFreq = 0; z_HPAmp = 0; z_FSlope = 0; z_ASlope = 0; z_SEntropy = 0; z_TEntropy = 0; z_STEntropy = 0; 
            for j = 1:sylnum
                source = [songsource '\' eval(['SylParam_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '(' num2str(index(j)) ').name'])];
                filenew = eval(['SylParam_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '(' num2str(index(j)) ').name']);
                filenew = regexprep(filenew, '.wav', ['_' SylNames(k) '.wav']);
                dest_app = [dest '\' SylNames(k)];
                if isdir(dest_app)
                else
                    mkdir(dest_app);
                end
                dest_app = [dest_app '\' filenew];
                copyfile(source, dest_app);
                eval(['FeatureSet_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '_' SylNames(k) '(j,:) = DataSet_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '(:,' num2str(index(j)) ');']);
                eval(['SylSet_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '_' SylNames(k) '(j,:) =  filenew;']);
                cd(dest);
                save(['SummaryData_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '_' SylNames(k) '.mat'], ['FeatureSet_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '_' SylNames(k)], ['SylSet_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '_' SylNames(k)]);
                %Copy feature set to temp arrays to calulate syllable drift
                z_Dur = z_Dur + eval(['SylParam_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '(' num2str(index(j)) ').z_Dur']);
                z_FFreq = z_FFreq + eval(['SylParam_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '(' num2str(index(j)) ').z_FFreq']);
                z_HPAmp = z_HPAmp + eval(['SylParam_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '(' num2str(index(j)) ').z_HPAmp']);
                z_FSlope = z_FSlope + eval(['SylParam_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '(' num2str(index(j)) ').z_FSlope']);
                z_ASlope = z_ASlope + eval(['SylParam_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '(' num2str(index(j)) ').z_ASlope']);
                z_SEntropy = z_SEntropy + eval(['SylParam_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '(' num2str(index(j)) ').z_SEntropy']);
                z_TEntropy = z_TEntropy + eval(['SylParam_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '(' num2str(index(j)) ').z_TEntropy']);
                z_STEntropy = z_STEntropy + eval(['SylParam_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '(' num2str(index(j)) ').z_STEntropy']);
            end
            %Substitute mean of z-scored syllable features for the next
            %day's target syllables.
            TestFeatures(1,k) = z_Dur/sylnum;
            TestFeatures(2,k) = z_FFreq/sylnum;
            TestFeatures(3,k) = z_HPAmp/sylnum;
            TestFeatures(4,k) = z_FSlope/sylnum;
            TestFeatures(5,k) = z_ASlope/sylnum;
            TestFeatures(6,k) = z_SEntropy/sylnum;
            TestFeatures(7,k) = z_TEntropy/sylnum;
            TestFeatures(8,k) = z_STEntropy/sylnum;
            %TestNames(k,:) = 'extracted'
            
        end
        eval(['SeedSyllables_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') ' =  TestFeatures;']);
        cd(dest);
        save(['SelectParam_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '.mat'], ['SeedSyllables_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd')]);
    end
    clear FeatureSet*
    clear DataSet*
    clear SylSet*
    date = datestr(addtodate(datenum(date), -1*interval, 'day'), 'yyyy-mm-dd');
    songsource = [parentsource '\' date];
    
    while exist([songsource '\' date '.mat'])~=2 & datenum(date, 'yyyy-mm-dd') >= datenum(earliest_date, 'yyyy-mm-dd')
        date = datestr(addtodate(datenum(date), -1, 'day'), 'yyyy-mm-dd');
        songsource = [parentsource '\' date];
    end
end
end