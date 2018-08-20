function VarSelectPCA(end_date, earliest_date, interval, test_syls)
%

%Ask user which file to find the dated folders of individuated syllables
%parentsource = uigetdir(pwd,'Source Folder');
parentsource = 'F:\Rd278 Stripped Syllables';

%Ask user where to store selected syllables
%parentdest = uigetdir(pwd,'Destination Folder');
parentdest = 'C:\Documents and Settings\Bird Brain\Desktop\Rd278 Sorted PCA 1.2';

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
    
TestNames = [];
TestFeatures = [];
for k=1:length(test_syls)
    TestNames(k,:) = eval(['SylParam_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '(' num2str(test_syls(k)) ').name']);
    TestFeatures(:,k) = eval(['DataSet_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd') '(:,' num2str(test_syls(k)) ')']);
end

while datenum(date, 'yyyy-mm-dd') >= datenum(earliest_date, 'yyyy-mm-dd')
    songsource = [parentsource '\' date];
    if isdir(songsource)
        load ([songsource '\' date '.mat']);
        dest = [parentdest '\' date];
        %The following crazy shit is required to get around MATLAB memory
        %constraints
        cur_DataSet = [TestFeatures eval(['DataSet_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd')])];
        [PCA_set, TestFeaturesPCA]=getPCA(cur_DataSet, size(TestFeatures,2));
        mem_index = [1:1000:length(PCA_set), length(PCA_set)];
        simfeat = [];
        for q = 1:(length(mem_index)-1)
            if q == (length(mem_index)-1)
                TempSet = [TestFeaturesPCA PCA_set(:,mem_index(q):(mem_index(q+1)))];
            else
                TempSet = [TestFeaturesPCA PCA_set(:,mem_index(q):(mem_index(q+1)-1))];
            end
            featdist = squareform(pdist(TempSet'));
            if q == 1
                simfeat = [simfeat featdist(1:size(TestFeaturesPCA,2),:)];
            else
                simfeat = [simfeat featdist(1:size(TestFeaturesPCA,2),size(TestFeaturesPCA,2)+1:end)];
            end
        end
        %jndDataSet = [TestFeatures eval(['DataSet_' datestr(datenum(date, 'yyyy-mm-dd'), 'yyyy_mm_dd')])];
        %simfeat=squareform(pdist(jndDataSet'));
        %***************
        dist_th = 1.2;   %This is the distance in feature space within which 2 syllables are "the same"; needs improvement.
        %***************
        if isdir(dest)
        else
            mkdir(dest);
        end
        for k=1:length(test_syls)
            index = find(0 < simfeat(k,length(test_syls)+1:end) & simfeat(k,length(test_syls)+1:end) <= dist_th);%-length(test_syls);
            testnum = min(length(index), 500);
            z_Dur = 0; z_FFreq = 0; z_HPAmp = 0; z_FSlope = 0; z_ASlope = 0; z_SEntropy = 0; z_TEntropy = 0; z_STEntropy = 0; 
            for j = 1:testnum
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
            TestFeatures(1,k) = z_Dur/testnum;
            TestFeatures(2,k) = z_FFreq/testnum;
            TestFeatures(3,k) = z_HPAmp/testnum;
            TestFeatures(4,k) = z_FSlope/testnum;
            TestFeatures(5,k) = z_ASlope/testnum;
            TestFeatures(6,k) = z_SEntropy/testnum;
            TestFeatures(7,k) = z_TEntropy/testnum;
            TestFeatures(8,k) = z_STEntropy/testnum;
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


