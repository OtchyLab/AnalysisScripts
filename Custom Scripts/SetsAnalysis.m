function [output] = SetsAnalysis
%Takes in the data created by Vernier GUI and generates summary data


%Get the Processed file to analyze
[temp,path] = uigetfile('C:\Users\Tim\Desktop\Ali Datasets Again\*.mat','Select the processed data file.','MultiSelect','on');
if size(temp,2)>1
    for i = 1:size(temp,2)
        fnames(i).name = char(temp{i});
    end
else
    fnames(1).name = temp;
end

%Predefine variables


for i = 1:length(fnames)
    %Load data
    load([path filesep fnames(i).name],'vernier')
    allVernier(i) = vernier;
    output(i).name = fnames(i).name;
    
    %Determine how many target intervals there are (i.e., syllable and gap)
    ints = vernier.targetEnd-vernier.targetStart;
    
    if ints == 1
        %Calculate audio stretch percent
        finalLength = vernier.audioPath(vernier.targetEnd,2)-vernier.audioPath(vernier.targetStart,2);
        startLength = vernier.audioPath(vernier.targetEnd,1)-vernier.audioPath(vernier.targetStart,1);
        output(i).audioStretch1 = 100*(finalLength-startLength)/startLength;
    
        %Calculate neuro stretch percent
        finalLength = vernier.neuroPath(vernier.targetEnd,2)-vernier.neuroPath(vernier.targetStart,2);
        startLength = vernier.neuroPath(vernier.targetEnd,1)-vernier.neuroPath(vernier.targetStart,1);
        output(i).neuroStretch1 = 100*(finalLength-startLength)/startLength;
        
        %Enter placeholder
        output(i).audioStretch2 = [];
        output(i).neuroStretch2 = [];
    elseif ints == 2
        %Calculate audio stretch percent for first interval
        finalLength = vernier.audioPath(vernier.targetStart+1,2)-vernier.audioPath(vernier.targetStart,2);
        startLength = vernier.audioPath(vernier.targetStart+1,1)-vernier.audioPath(vernier.targetStart,1);
        output(i).audioStretch1 = 100*(finalLength-startLength)/startLength;
    
        %Calculate neuro stretch percent for first interval
        finalLength = vernier.neuroPath(vernier.targetStart+1,2)-vernier.neuroPath(vernier.targetStart,2);
        startLength = vernier.neuroPath(vernier.targetStart+1,1)-vernier.neuroPath(vernier.targetStart,1);
        output(i).neuroStretch1 = 100*(finalLength-startLength)/startLength;
        
        %Calculate audio stretch percent for first interval
        finalLength = vernier.audioPath(vernier.targetEnd,2)-vernier.audioPath(vernier.targetStart+1,2);
        startLength = vernier.audioPath(vernier.targetEnd,1)-vernier.audioPath(vernier.targetStart+1,1);
        output(i).audioStretch2 = 100*(finalLength-startLength)/startLength;
    
        %Calculate neuro stretch percent for first interval
        finalLength = vernier.neuroPath(vernier.targetEnd,2)-vernier.neuroPath(vernier.targetStart+1,2);
        startLength = vernier.neuroPath(vernier.targetEnd,1)-vernier.neuroPath(vernier.targetStart+1,1);
        output(i).neuroStretch2 = 100*(finalLength-startLength)/startLength;
    else
        print('too many intervals')
        return
    end
    
    %Set up the intervals of interest
    offset = vernier.audioPath(vernier.targetStart,1)+vernier.timeVect(1);
    targetOnlyR = (vernier.audioPath(vernier.targetStart,1)-offset):(vernier.audioPath(vernier.targetEnd,1)-offset);
    nontargetR = (vernier.audioPath(2)-offset):(vernier.audioPath(vernier.targetStart,1)-offset);
    if (vernier.audioPath(vernier.targetEnd,1)-offset+100) > length(vernier.unwarpCorr)
        targetPlusR =(vernier.audioPath(vernier.targetStart,1)-offset):length(vernier.unwarpCorr);
    else    
        targetPlusR =(vernier.audioPath(vernier.targetStart,1)-offset):(vernier.audioPath(vernier.targetEnd,1)-offset+100);
    end
    
    %Calculate correlation for target only
    output(i).targetOnlyUnCorr = mean(vernier.unwarpCorr(targetOnlyR));
    output(i).targetOnlyWaCorr = mean(vernier.warpCorr(targetOnlyR));
    
    %Calculate correlation for target plus 100
    output(i).targetPlusUnCorr = mean(vernier.unwarpCorr(targetPlusR));
    output(i).targetPlusWaCorr = mean(vernier.warpCorr(targetPlusR));
    
    %Calculate correlation for non-target
    output(i).nontargetUnCorr = mean(vernier.unwarpCorr(nontargetR));
    
    clear('vernier')
end


%Loop through the output structure to collect the data of interest
summary.audioStretch = [];
summary.neuroStretch = [];
for i = 1:length(fnames)
    summary.un_targetOnly(i) = output(i).targetOnlyUnCorr;
    summary.un_targetPlus(i) = output(i).targetPlusUnCorr;
    
    summary.warp_targetOnly(i) = output(i).targetOnlyWaCorr;
    summary.warp_targetPlus(i) = output(i).targetPlusWaCorr;
    
    summary.un_nonTarget(i) = output(i).nontargetUnCorr;
    
    summary.audioStretch = [summary.audioStretch;output(i).audioStretch1];
    summary.neuroStretch = [summary.neuroStretch;output(i).neuroStretch1];
    
    if ~isempty(output(i).audioStretch2) || ~isempty(output(i).neuroStretch2)
        summary.audioStretch = [summary.audioStretch;output(i).audioStretch2];
        summary.neuroStretch = [summary.neuroStretch;output(i).neuroStretch2];
    end
end

%Calculate sample stats
summary.un_targetOnlyM = mean(summary.un_targetOnly);
summary.un_targetOnlyS = std(summary.un_targetOnly)/sqrt(i);

summary.warp_targetOnlyM = mean(summary.warp_targetOnly);
summary.warp_targetOnlyS = std(summary.warp_targetOnly)/sqrt(i);

summary.un_targetPlusM = mean(summary.un_targetPlus);
summary.un_targetPlusS = std(summary.un_targetPlus)/sqrt(i);

summary.warp_targetPlusM = mean(summary.warp_targetPlus);
summary.warp_targetPlusS = std(summary.warp_targetPlus)/sqrt(i);

summary.un_nonTargetM = mean(summary.un_nonTarget);
summary.un_nonTargetS = std(summary.un_nonTarget)/sqrt(i);

summary.stretch_R = corr(summary.audioStretch,summary.neuroStretch);

%Calculate p-values
[h,summary.pVal_targetOnly2nontarget] = ttest2(summary.un_targetOnly,summary.un_nonTarget) %for pCAF
[h,summary.pVal_targetPlus2nontarget] = ttest2(summary.un_targetPlus,summary.un_nonTarget) %for tCAF
[h,summary.pVal_warptargetPlus2nontarget] = ttest2(summary.warp_targetPlus,summary.un_nonTarget) %for tCAF

%Plot the mean and sem bars for summary Stats for tCAF figure
figure
hold on
bar(0,summary.un_targetPlusM,'FaceColor','none','EdgeColor',[0 0 0],'BarWidth',0.6);
bar(2,summary.warp_targetPlusM,'FaceColor','none','EdgeColor',[0 0 0],'BarWidth',0.6);
bar(1,summary.un_nonTargetM,'FaceColor','none','EdgeColor',[0 0 0],'BarWidth',0.6);
errorbar(0,summary.un_targetPlusM,summary.un_targetPlusS,'Marker','.','LineStyle','none','Color',[0 0 0])
errorbar(2,summary.warp_targetPlusM,summary.warp_targetPlusS,'Marker','.','LineStyle','none','Color',[0 0 0])
errorbar(1,summary.un_nonTargetM,summary.un_nonTargetS,'Marker','.','LineStyle','none','Color',[0 0 0])
axis tight
xlim([-.5,2.5])
ylim([0,1])
set(gca,'XTick', [],'YTick',[0,.5,1])
set(gca,'TickDir','out','Box','off')
ylabel('Mean Correlation')
title('Correlation in HVC Activity Pre and Post tCaf')

%Plot the mean and sem bars for summary Stats for pCAF figure
figure
hold on
bar(0,summary.un_targetOnlyM,'FaceColor','none','EdgeColor',[0 0 0],'BarWidth',0.6);
bar(1,summary.un_nonTargetM,'FaceColor','none','EdgeColor',[0 0 0],'BarWidth',0.6);
errorbar(0,summary.un_targetOnlyM,summary.un_targetOnlyS,'Marker','.','LineStyle','none','Color',[0 0 0])
errorbar(1,summary.un_nonTargetM,summary.un_nonTargetS,'Marker','.','LineStyle','none','Color',[0 0 0])
axis tight
xlim([-.5,2.5])
ylim([0,1])
set(gca,'XTick', [],'YTick',[0,.5,1])
set(gca,'TickDir','out','Box','off')
ylabel('Mean Correlation')
title('Correlation in HVC Activity Pre and Post pCaf')

%Plot the scatter of the strecth intervals
figure
hold on
line([-50,50],[-50,50],'Color','k','LineStyle','--','LineWidth',1.5)
scatter(summary.audioStretch,summary.neuroStretch,'or')
xlabel('Song Warping (%)')
ylabel('Neural Warping (%)')

a=1;
%Get location to save output data to file.
% [filename,pathLoc] = uiputfile('C:\Users\Tim\Desktop\Ali Datasets Again\*.mat');
% save([pathLoc filename],'allVernier','output','summary')

