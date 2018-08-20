function xcorrMatrixBATCH()

%Get the Processed file to analyze
[temp,path] = uigetfile('*.mat','Select the processed data file.','MultiSelect','on');

if size(temp,2)>1
    for i = 1:size(temp,2)
        fnames(i).name = char(temp{i});
    end
else
    fnames.name = temp;
end

for i = 1:length(fnames)
    load([path filesep fnames(i).name],'data')
    nPowerMat = data.neuroPower;
    clear('data')

    %Generate covariance matrices
    [neuroCov,neuroP] = corrcoef(nPowerMat');
    neuroCovArray = tril(neuroCov,-1);
    neuroCovArray = neuroCovArray(neuroCovArray~=0);
    name_neuroCov{i} = fnames(i).name;
    m_neuroCov(i) = mean(neuroCovArray);
    std_neuroCov(i) = std(neuroCovArray);
    
    t = regexp(fnames(i).name,' ','split');
    time = t{4}(1:end-6);
    if ~strcmp(time(1),'n')
        lag_neuroCov(i) = str2num(time);
    else
        lag_neuroCov(i) = -1*str2num(time(4:end));
    end
    
    %Show data on axes
    fig1 = figure;
    subplot(2,1,1)
    h=imagesc(neuroCov);
    title(['Correlation Matrix for ' fnames(i).name])

    subplot(2,1,2)
    ecdf(neuroCovArray);

    %Save image and data to file
    save([path 'Data\' fnames(i).name(1:end-4), ' data.mat'],'neuroCov','neuroCovArray')
    saveas(fig1,[path 'Images\' fnames(i).name(1:end-4), '.fig'])
    close(fig1)
    clear('neuroCov','neuroP','neuroCovArray','nPowerMat')
end

%Resort all of the data by lags
[lag_neuroCov, indx] = sort(lag_neuroCov);
name_neuroCov = name_neuroCov(indx);
m_neuroCov = m_neuroCov(indx);
std_neuroCov = std_neuroCov(indx);

%Save the summary to file
save([path 'Data\' fnames(i).name(1:17), ' summary.mat'],'name_neuroCov','lag_neuroCov','m_neuroCov','std_neuroCov')



