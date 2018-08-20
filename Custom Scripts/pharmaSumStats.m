function pharmaSumStats
%4/7/15
%Updated from the original version (inactSumStats.m) to remove interpolation of heat maps, reference
%the original experiment type. Principally for the PHARMACOLOGICAL Nif lesions.
%Takes the output of the lesionAnalysisV2.m scripts (across several days, and all in the same folder) 
%and runs summary stats on the different conditions. All output data is saved to an on-disk mat-file that include
%bird names, dates, conditions, and results.

%Set the output file location
outputLocation = 'C:\Users\Tim\Desktop\Nif Project Figures\Pharma\';
outName = 'Pharma Summary Stats EMD 150513.mat';

%Select the dataset to use
%[dataName, dataFolder] = uigetfile('C:\Users\Tim\Desktop\Nif Project Figures\Inactivations\Grn114\*.mat');
% dataFolder = 'C:\Users\Tim\Desktop\Nif Project Figures\Pharma\Grn089\';
% dataName = 'Grn089_140625_pharmaDataset.mat';
% load([dataFolder, filesep, dataName]);
% close all; %close the windows

%Ask user for data directory
dataFolder = uigetdir('C:\Users\Tim\Desktop\Nif Project Figures\Pharma\');

%Get all data files that match the criterion
crit = '*pharmaDataset.mat';
filelist = dir([dataFolder, filesep, crit]);

%Cycle through the filelist and sequentially process the datasets from the whole folder
data = [];
output = [];
for i = 1:length(filelist)
    %Select the file name
    dataName = filelist(i).name;

    %Load the data file
    t = load([dataFolder, filesep, dataName]);
    
    %Parse name to grab bird name, date and type of experiment
    sp = regexp(dataName,'_', 'split');
    output.name = char(sp{1});
    
    %Correct for L-5d to L-4d switch so rest of code works
    if strcmp(char(sp{1}), 'Grn091') && t.conIdx==2
        t.conIdx = 1;
        t.offset = -5;
        display('Made the Grn091 switch!')
    end
    
    output.date{t.conIdx} = char(sp{2});
    output.filename{t.conIdx} = [dataFolder, filesep, dataName];
    output.type{t.conIdx} = t.conditions{t.conIdx};
    output.offset{t.conIdx} = t.offset;
    output.conditions = t.conditions;
    
    % Capture distributions for later stats
    output.pdfY(t.conIdx,:) = t.pdfY;
    output.N(t.conIdx,:,:) = t.N;
%     output.Ncnts(t.conIdx,:,:) = t.Ncnts;
    
    %Filter the original PDFs with a gaussian and store in original variable for ease
    output.filter.size = 7;
    output.filter.sigma = 3;
    output.filter.h = fspecial('gaussian', output.filter.size, output.filter.sigma);
    output.Nsm(t.conIdx,:,:) = imfilter(squeeze(output.N(t.conIdx,:,:)), output.filter.h, 'replicate');
%     figure; imagesc(squeeze(output.Nsm(t.conIdx,:,:))); colormap(jet);
    close all; %close the windows
    
    %Clear the dataset
    clear('t')
    
end

%%%%%%%%%%%%%%%
% Averaged Pre
%%%%%%%%%%%%%%%
%Compare day of lesion data (col 1)  to all others (col 2)
cmpMat = [6*ones(11,1),[1:11]'];

%Copy to output structure
output.cmpMat = cmpMat;

%pre-load costmatrix
load('CostMatrix 281x161.mat');

%Loop through the comparisons and calculate stats
for i = 1:size(cmpMat,1)
    
    if isempty(output.offset{cmpMat(i, 2)})
        %If this subset is empty, do nothing???
        
    else
        %Otherwise, do the processing on thie dataset.
        
        %%%%%%%%%%%%%%%%%%%%
        % 1-D measures
        %%%%%%%%%%%%%%%%%%%%

        %Simplify variables for reading
        %To define the pre template, take the mean of up to 2 non-empty days pre-lesion
%         preY = [];
%         for j = 4:6
%             if ~isempty(output.offset{j})
%                 preY = [preY; output.pdfY(j,:)];
%             end
%         end
%         pdfY1 = mean(preY,1);
%         output.preY = pdfY1;
%         
%         %Compare to
%         pdfY2 = output.pdfY(cmpMat(i, 2),:);
% 
%         %Simple Pearson Linear Correlation 
%         output.D1corr(i) = corr(pdfY1', pdfY2');
% 
%         % Euclidean Distance in 1D space
%         output.D1dist(i) = pdist([pdfY1; pdfY2], 'euclidean');
% 
%         % Cosine (i.e., angle between the normalized pdf vectors)
%         output.D1cos(i) = acos(dot(pdfY1,pdfY2)/(norm(pdfY1)*norm(pdfY2)));
% 
%         % Bhattacharyya distance
%         output.D1bhatt(i) = -log(sum(sqrt(pdfY1.*pdfY2)));
% 
%         % Hellinger distance
%         output.D1hell(i) = sqrt(2*(sum((sqrt(pdfY1) - sqrt(pdfY2)).^2)));
% 
%         %Replace zero values with very small values to avoid singlularity points (i.e., 0*log(0) = NaN)
%         zeta = 10^-100;
%         pdfY1zeta = pdfY1; pdfY1zeta(pdfY1==0) = zeta;
%         pdfY2zeta = pdfY2; pdfY2zeta(pdfY2==0) = zeta;
% 
%         % Kullback-Leibler divergence (relative entropy)
%         output.D1KL(i) = 0;
% 
%         %Shannon entropy of the distribution (no comparison between treatments)
%          output.D1se(i) = sum(pdfY1zeta.*log(pdfY2zeta));

        %%%%%%%%%%%%%%%%%%%%
        % 2-D measures
        %%%%%%%%%%%%%%%%%%%%
        %Simplify variables for reading
        %To define the pre template, take the mean of up to 3 non-empty days pre-lesion
        preN = []; preNcnts = [];
        for j = 4:5
            if ~isempty(output.offset{j})
                preN = [preN; output.Nsm(j,:,:)]; %Use filtered version
%                 preNcnts = [preNcnts; output.Ncnts(j,:,:)];
            end
        end
        N1 = squeeze(mean(preN,1));
        output.preN = N1;
        
%         N1cnts = squeeze(sum(preNcnts,1));
%         output.preNcnts = N1cnts;
        
           %Simplify variables for reading
        N2 = squeeze(output.Nsm(cmpMat(i, 2),:,:)); %Use filtered version
%         N2cnts = squeeze(output.Ncnts(cmpMat(i, 2),:,:));

%         %Simple Pearson Linear Correlation 
%         output.D2corr(i) = corr(N1(:), N2(:));
% 
%         % Euclidean Distance in 1D space
%         output.D2dist(i) = pdist([N1(:)'; N2(:)'], 'euclidean');

        % Cosine (i.e., angle between the normalized pdf vectors)
%         output.D2cos(i) = acos(dot(N1(:),N2(:))/(norm(N1(:))*norm(N2(:))));
        output.D2emd(i) = calcEMD(N1,N2,D); %Temp substitute the EMD

%         % Bhattacharyya distance
%         output.D2bhatt(i) = -log(sum(sqrt(N1(:).*N2(:))));
% 
%         % Hellinger distance
%         output.D2hell(i) = sqrt(2*(sum((sqrt(N1(:)) - sqrt(N2(:))).^2)));
% 
%         %Replace zero values with very small values to avoid singlularity points (i.e., 0*log(0) = NaN)
%         zeta = 10^-100;
%         N1zeta = N1(:); N1zeta(N1(:)==0) = zeta;
%         N2zeta = N2(:); N2zeta(N2(:)==0) = zeta;
% 
%         % Kullback-Leibler divergence (relative entropy)
% %         output.D2KL(i) = sum(N1zeta.*log(N1zeta./N2zeta));
%         output.D2KL(i) = getKL(N1cnts,N2cnts);
% 
%         %Shannon entropy of the distribution (no comparison between treatments)
%         output.D2se(i) = sum(N1zeta.*log(N2zeta));
    end
end

%%%%%%%%%%%%%%%%%%%%%
% Save to output
%%%%%%%%%%%%%%%%%%%%%
m = exist([outputLocation, outName]);
if m ==2 
    %File already exists
    clear('pharmaStats');
    load([outputLocation, outName], 'pharmaStats')
    pharmaStats(end+1) = output;
else
    %No file yet created
    pharmaStats = output;
end

%Save the updated data to file
save([outputLocation, outName], 'pharmaStats')

display(['done with ' output.name])

clear all



