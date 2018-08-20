function inactSumStats
%Takes the output of the inactAnalysis scripts (the complete set of working data) and runs summary stats on the different
%conditions. All output data is saved to an on-disk mat-file that include bird names, dates, conditions, and results.

clear all; close all;

%Set the output file location
outputLocation = 'C:\Users\Tim\Desktop\';
outName = 'Inactivation Summary Stats.mat';

%Select the dataset to use
[dataName, dataFolder] = uigetfile('C:\Users\Tim\Desktop\Nif Project Figures\Inactivations\Grn114\*.mat');
% dataFolder = 'C:\Users\Tim\Desktop\Nif Project Figures\Inactivations\Grn108\';
% dataName = 'Grn108_140828_inactDataset.mat';
load([dataFolder, dataName]);
close all; %close the windows

% injType = 'Inact';
injType = 'PBS';
% injType = 'Elev';

%Parse name to grab bird name, date and type of experiment
sp = regexp(dataName,'_', 'split');
output.name = char(sp{1});
output.date = char(sp{2});
output.filename = [dataFolder, dataName];
output.type = injType;

%Comparison labels
output.labels = [{'pre-inj'}, {'inj-post'}, {'pre-post'}];

%%%%%%%%%%%%%%%%%%%%
% 1-D measures
%%%%%%%%%%%%%%%%%%%%
% 
% %Simple Pearson Linear Correlation 
% output.D1corr = [corr(preY', injY'), corr(injY', postY'), corr(preY', postY')];
% 
% % Euclidean Distance in 1D space
% d = pdist([preY; injY; postY],'euclidean');
% output.D1dist = [d(1), d(3), d(2)]; %rearange to keep consistent
% 

%%%%%%%%%%%%%%%%%%%%
% 2-D measures
%%%%%%%%%%%%%%%%%%%%

%Interpolate and normalizes all the heat maps
[edgesX, edgesY, N, ~] = ndhist(sylFeat(preInx), sylDur(preInx), 'bins', 6, 'filter', 'axis', [-3.5 0 0 400]);
interpN = interp2(edgesX, edgesY, N, linspace(-4, 0,100)', linspace(0, 400, 150));
interpN(isnan(interpN)) = 0;
norm_interpPreN = interpN./(sum(interpN(:)));

[edgesX, edgesY, N, ~] = ndhist(sylFeat(injInx), sylDur(injInx), 'bins', 6, 'filter', 'axis', [-3.5 0 0 400]);
interpN = interp2(edgesX, edgesY, N, linspace(-4, 0,100)', linspace(0, 400, 150));
interpN(isnan(interpN)) = 0;
norm_interpInjN = interpN./(sum(interpN(:)));

[edgesX, edgesY, N, ~] = ndhist(sylFeat(postInx), sylDur(postInx), 'bins', 6, 'filter', 'axis', [-3.5 0 0 400]);
interpN = interp2(edgesX, edgesY, N, linspace(-4, 0,100)', linspace(0, 400, 150));
interpN(isnan(interpN)) = 0;
norm_interpPostN = interpN./(sum(interpN(:)));

close all; %close the windows

%Simple Pearson Linear Correlation 
output.D2corr = [corr(norm_interpPreN(:), norm_interpInjN(:)), corr(norm_interpInjN(:), norm_interpPostN(:)), corr(norm_interpPreN(:), norm_interpPostN(:))];

%Euclidean Distance in 2D space
d = pdist([norm_interpPreN(:)'; norm_interpInjN(:)'; norm_interpPostN(:)'], 'euclidean');
output.D2dist = [d(1), d(3), d(2)]; %rearange to keep consistent

%%%%%%%%%%%%%%%%%%%%%
% Save to output
%%%%%%%%%%%%%%%%%%%%%

m = exist([outputLocation, outName]);
if m ==2 
    %File already exists
    clear('inactStats');
    load([outputLocation, outName], 'inactStats')
    inactStats(end+1) = output;
else
    %No file yet created
    inactStats = output;
end

%Save the updated data to file
save([outputLocation, outName], 'inactStats')

clear all

display('done')

