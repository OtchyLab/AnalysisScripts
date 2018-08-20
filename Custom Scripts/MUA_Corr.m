function MUA_Corr
%The script takes the audio and neural alignments (created by StretchEm)
%and generates the pairwise correlation matrices for both.  This is
%substantially the same as the XCorr function within StretchEm.
clear all

loc = 'C:\Users\Tim\Desktop\DirUndir XCorr Sets\Alignment Sets';
cd(loc)
files = dir('P*.mat');
numFiles = length(files);

maxLen = 0;
pntr = 0;
for i = 1:numFiles
    %Load the structures I care about
    load(files(i).name,'data')
    
    corrVar{i}.name = files(i).name;
    
    %Find the break point for the dir/und
    RemapIndx = data.directRemapInx;
    breaks = find(diff(RemapIndx)<0,1,'first')+1;
    corrVar{i}.breaks = breaks;
    
    if i == numFiles %Special case for this dataset!!!
        breaks = 31;
        corrVar{i}.breaks = breaks;
    end
    
    sc_audio = data.audioPower(RemapIndx,50:end-50);
    sc_neuro = data.neuroPower(RemapIndx,(.05*44150):end-(.05*44150));

    %Generate covariance matrices for each
    [audCov,audP] = corrcoef(sc_audio');
    [neuroCov,neuroP] = corrcoef(sc_neuro');
    
    %audCov(audCov==1) = NaN;
    %neuroCov(neuroCov==1) = NaN;

    %Show data on axes
    h = figure;
    subplot(2,1,1)
    corrVar{i}.audioMat = audCov./(mean(mean(audCov)));
    h1=imagesc(corrVar{i}.audioMat);
    hold on
    line([0 length(RemapIndx)],[breaks breaks],'Color','w','LineWidth',2)
    line([breaks breaks],[0 length(RemapIndx)],'Color','w','LineWidth',2)
    hold off
    title('Rendition-by-rendition cross correlation of Song Power')
    set(gca,'TickDir','out','Box','off')

    subplot(2,1,2)
    corrVar{i}.neuroMat = neuroCov./(mean(mean(neuroCov)));
    h2=imagesc(corrVar{i}.neuroMat);
    hold on
    line([0 length(RemapIndx)],[breaks breaks],'Color','w','LineWidth',2)
    line([breaks breaks],[0 length(RemapIndx)],'Color','w','LineWidth',2)
    hold off
    title('Rendition-by-rendition cross correlation of HVC Power')
    set(gca,'TickDir','out','Box','off')
    
    %Calculate mean and std for the directed and undirected xcorrs
    nUndir = neuroCov(1:breaks-1,1:breaks-1);
    nUndir = nUndir(:);
    corrVar{i}.nUndir = nUndir(nUndir~=1);

    nDir = neuroCov(breaks:end,breaks:end);
    nDir = nDir(:);
    corrVar{i}.nDir = nDir(nDir~=1);

    %Calculate Stats on the datasets
    corrVar{i}.neuro_undir_m = mean(nUndir);
    corrVar{i}.neuro_undir_std = std(nUndir);

    corrVar{i}.neuro_dir_m = mean(nDir);
    corrVar{i}.neuro_dir_std = std(nDir);
    
    %Generate eCDFs for the two distributions
    m = figure;
    hold on
    h3 = cdfplot(corrVar{i}.nUndir);
    h4 = cdfplot(corrVar{i}.nDir);
    
    [r,corrVar{i}.ttestp]=ttest2(corrVar{i}.nUndir,corrVar{i}.nDir);
    text(0.15,.65,['p<' num2str(corrVar{i}.ttestp,2)],'EdgeColor','k','FontSize',20)
    set(h3,'DisplayName','Undir','LineWidth',2,'Color','b')
    set(h4,'DisplayName','Dir','LineWidth',2,'Color','r')
    set(gca,'TickDir','out','Box','off')
    grid(gca,'off');
    xlim([0,1]);xlabel('Xcorr-Value')
    ylim([0,1]);ylabel('P(x)')
    title('CDF for HVC Power in Undir and Dir')

    %Create summary figure with means and std
    fig666 = figure(666);
    hold on
    f = errorbar([1,3],...
             [corrVar{i}.neuro_undir_m,corrVar{i}.neuro_dir_m],...
             [corrVar{i}.neuro_undir_std,corrVar{i}.neuro_dir_std],...
             ':.','LineWidth',1,'DisplayName',corrVar{i}.name(1:end-4));
         
             %Create summary figure with means and std
    fig667 = figure(667);
    hold on
    f = errorbar([1,3],...
             [corrVar{i}.neuro_undir_m/corrVar{i}.neuro_undir_m,corrVar{i}.neuro_dir_m/corrVar{i}.neuro_undir_m],...
             [corrVar{i}.neuro_undir_std/corrVar{i}.neuro_undir_m,corrVar{i}.neuro_dir_std/corrVar{i}.neuro_undir_m],...
             ':o','LineWidth',1,'DisplayName',corrVar{i}.name(1:end-4));
        
    %Save plots to file
    clear('data')
    saveas(h,[loc '\XCorr\' corrVar{i}.name(1:14), corrVar{i}.name(22:25) ' xcorr.fig'])
    saveas(h,[loc '\XCorr\' corrVar{i}.name(1:14), corrVar{i}.name(22:25) ' xcorr.tif'])
    saveas(m,[loc '\XCorr\' corrVar{i}.name(1:14), corrVar{i}.name(22:25) ' CDF.fig'])
    saveas(m,[loc '\XCorr\' corrVar{i}.name(1:14), corrVar{i}.name(22:25) ' CDF.tif'])
    close(h)
    close(m)
    
end

%Statistics across days
b1 = [1,2,3]; %Pur639
b2 = [4,5]; %Pur683
b3 = [6,7,8]; %Pur692
b4 = [9,10]; %Pur755

neuro_undir_m = getAnnotationVector(corrVar,'neuro_undir_m');
%neuro_undir_std = getAnnotationVector(corrVar,'neuro_undir_std');
neuro_dir_m = getAnnotationVector(corrVar,'neuro_dir_m');
%neuro_dir_std = getAnnotationVector(corrVar,'neuro_dir_std');

figure(666)
hold on
f = errorbar([1,3],...
             [mean(neuro_undir_m(b1)),mean(neuro_dir_m(b1))],...
             [std(neuro_undir_m(b1)),std(neuro_dir_m(b1))],...
             '-sb','LineWidth',2,'DisplayName',corrVar{1}.name(1:end-10));
f = errorbar([1,3],...
             [mean(neuro_undir_m(b2)),mean(neuro_dir_m(b2))],...
             [std(neuro_undir_m(b2)),std(neuro_dir_m(b2))],...
             '-sr','LineWidth',2,'DisplayName',corrVar{4}.name(1:end-10));
f = errorbar([1,3],...
             [mean(neuro_undir_m(b3)),mean(neuro_dir_m(b3))],...
             [std(neuro_undir_m(b3)),std(neuro_dir_m(b3))],...
             '-sg','LineWidth',2,'DisplayName',corrVar{6}.name(1:end-10));
f = errorbar([1,3],...
             [mean(neuro_undir_m(b4)),mean(neuro_dir_m(b4))],...
             [std(neuro_undir_m(b4)),std(neuro_dir_m(b4))],...
             '-sm','LineWidth',2,'DisplayName',corrVar{9}.name(1:end-10));

%Format Summary figure
title('Mean +/- Std of HVC Power Pairwise Correlations')
ylabel('Mean Pairwise Correlation')
set(gca,'XTick',[1,3],'XTickLabel',{char('UNDIR'),char('DIR')})
set(gca,'TickDir','out','Box','off')

figure(667)
hold on
f = errorbar([1,3],...
             [mean(neuro_undir_m(b1))/mean(neuro_undir_m(b1)),mean(neuro_dir_m(b1))/mean(neuro_undir_m(b1))],...
             [std(neuro_undir_m(b1))/mean(neuro_undir_m(b1)),std(neuro_dir_m(b1))/mean(neuro_undir_m(b1))],...
             '-sb','LineWidth',2,'DisplayName',corrVar{1}.name(1:end-10));
f = errorbar([1,3],...
             [mean(neuro_undir_m(b2))/mean(neuro_undir_m(b2)),mean(neuro_dir_m(b2))/mean(neuro_undir_m(b2))],...
             [std(neuro_undir_m(b2))/mean(neuro_undir_m(b2)),std(neuro_dir_m(b2))/mean(neuro_undir_m(b2))],...
             '-sr','LineWidth',2,'DisplayName',corrVar{4}.name(1:end-10));
f = errorbar([1,3],...
             [mean(neuro_undir_m(b3))/mean(neuro_undir_m(b3)),mean(neuro_dir_m(b3))/mean(neuro_undir_m(b3))],...
             [std(neuro_undir_m(b3))/mean(neuro_undir_m(b3)),std(neuro_dir_m(b3))/mean(neuro_undir_m(b3))],...
             '-sg','LineWidth',2,'DisplayName',corrVar{6}.name(1:end-10));
f = errorbar([1,3],...
             [mean(neuro_undir_m(b4))/mean(neuro_undir_m(b4)),mean(neuro_dir_m(b4))/mean(neuro_undir_m(b4))],...
             [std(neuro_undir_m(b4))/mean(neuro_undir_m(b4)),std(neuro_dir_m(b4))/mean(neuro_undir_m(b4))],...
             '-sm','LineWidth',2,'DisplayName',corrVar{9}.name(1:end-10));

%Format Summary figure
title('Mean +/- Std of HVC Power Pairwise Correlations')
ylabel('Mean Pairwise Correlation')
set(gca,'XTick',[1,3],'XTickLabel',{char('UNDIR'),char('DIR')})
set(gca,'TickDir','out','Box','off')

%Save data to file
savename = [loc '\XCorr\' 'XcorrVars.mat'];
save(savename,'corrVar')

saveas(fig666,[loc '\XCorr\' 'XcorrSummary.fig'])
saveas(fig666,[loc '\XCorr\' 'XcorrSummary.tif'])
saveas(fig667,[loc '\XCorr\' 'NXcorrSummary.fig'])
saveas(fig667,[loc '\XCorr\' 'NXcorrSummary.tif'])