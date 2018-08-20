function [neuroCov,neuroCovArray] = xcorrMatrix()

%Get the Processed file to analyze
[file,path] = uigetfile('*.mat','Select the processed data file.');
load([path filesep file],'data')
nPowerMat = data.neuroPower;
clear('data')

%Generate covariance matrices
[neuroCov,neuroP] = corrcoef(nPowerMat');
neuroCovArray = tril(neuroCov,-1);
neuroCovArray = neuroCovArray(neuroCovArray~=0);

%Show data on axes
fig1 = figure;
subplot(2,1,1)
h=imagesc(neuroCov);
title(['Correlation Matrix for ' file])

subplot(2,1,2)
ecdf(neuroCovArray);

%Save image and data to file
save([path 'Data\' file(1:end-4), 'data.mat'],'neuroCov','neuroCovArray')
saveas(fig1,[path 'Images\' file(1:end-4), '.fig'])
