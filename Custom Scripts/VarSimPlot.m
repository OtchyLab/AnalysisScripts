function VarSimPlot(end_date, earliest_date, interval, DOBstr)
%Designed to operate on the output/saved files from the VarSelect routine.
%Function assumes a parent folder with sorted syllables arranged in dated
%files (yyyy-mm-dd).  In particular, .mat files named
%SummaryData_date_syltyp.mat and SelectParam_date.mat with Syllable
%information.

%Ask user which file to find the dated folders of individuated syllables
%parentsource = uigetdir(pwd,'Source Folder');
parentsource = 'C:\Documents and Settings\Bird Brain\Desktop\Rd278 Sorted 1.3';

date = end_date;
DOB=datenum(DOBstr, 'yyyy-mm-dd');

SylNames = ['A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'; 'K'; 'L'; 'M'; 'O'; 'P'];
SylSimAll = [];
SylSimAll2 = [];
SylAveAll = [];
tempSylList = [];
tempSylAve = [];

songsource = [parentsource '\' date];
underdate=regexprep(date, '-', '_');

while exist([songsource '\SummaryData_' underdate '_A.mat'])~=2
    date = datestr(addtodate(datenum(date), -1, 'day'), 'yyyy-mm-dd');
    underdate=regexprep(date, '-', '_');
    songsource = [parentsource '\' date];
end

while datenum(date, 'yyyy-mm-dd') >= datenum(earliest_date, 'yyyy-mm-dd')
    songsource = [parentsource '\' date];
    underdate=regexprep(date, '-', '_');
    curDateNum=datenum(date, 'yyyy-mm-dd');
    DPH = curDateNum-DOB;
    tempSylAve(1) = DPH;
    cd(songsource);
    files = dir('SummaryData*.mat');
    for k=1:length(files)
        load ([songsource '\' files(k).name], 'Feature*'); %Load all of the .mat files corresponding to the sorted Syls
        sylDist = pdist(eval(['FeatureSet_' underdate '_' SylNames(k)]));
        aveDist = mean(sylDist);
        numSyls = length(sylDist);
        tempSylAve(k+1) = aveDist;
        tempSylList(1:numSyls,1) = DPH; tempSylList(1:numSyls,k+1) = sylDist'; 
    end
    clear FeatureSet*

    if length(SylSimAll)>3000000
        SylSimAll2 = [SylSimAll2; tempSylList];
    else
        SylSimAll = [SylSimAll; tempSylList];
    end
    SylAveAll = [SylAveAll; tempSylAve, mean(tempSylAve(2:end))];

    date = datestr(addtodate(datenum(date), -1*interval, 'day'), 'yyyy-mm-dd');
    underdate=regexprep(date, '-', '_');
    songsource = [parentsource '\' date];
    while exist([songsource '\SummaryData_' underdate '_A.mat'])~=2
        date = datestr(addtodate(datenum(date), -1, 'day'), 'yyyy-mm-dd');
        underdate=regexprep(date, '-', '_');
        songsource = [parentsource '\' date];
    end
end

%Plot similarity/variability data for each syllable on an axis (sim vs
%DPH
syltypes = size(SylAveAll,2)-2;

%**************
zero_sim = 3.75; %Threshhold for 0% similarity--derived from the fairly 
                 %consistent mean distance between all syllables for a day
                 %Experimentally determined value
%**************

    for k=1:syltypes
        subplot(syltypes,2,k*2-1);
        sim1 = (zero_sim-SylSimAll(:,k+1)).*100./zero_sim; %Scales similarity to a 0-100% similarity scalar.
        sim2 = (zero_sim-SylSimAll2(:,k+1)).*100./zero_sim; %Scales similarity to a 0-100% similarity scalar.
        plot(SylSimAll(:,1), sim1, '.');
        hold on
        plot(SylSimAll2(:,1), sim2, '.');
        xlim([35 max(SylSimAll(1,1))+15]);
        ylim([0 100]);
        %ylabel('Similarity/Variability (ADU)');
        %xlabel('Days Post Hatch');
        hold on;

        %Plot variability on left-hand axes
        subplot(syltypes,2,k*2-1);
        sim = (zero_sim-SylAveAll(:,k+1)).*100./zero_sim;
        var = (100-sim).*100./50; %Calcs variability from similarity according to (Sim{max}-sim{i})/(Sim{max}-Sim{min})
        plot(SylAveAll(:,1), var, 'r');
        hold on

        %Plot variability (lightly)on right hand axis
        axes_sel=[1:syltypes]*2;
        subplot(syltypes,2,axes_sel);
        plot(SylAveAll(:,1), var, 'g');
        xlim([35 max(SylSimAll(1,1))+15]);
        ylim([0 100]);
        %ylabel('Similarity (ADU)');
        %xlabel('Days Post Hatch');
        hold on;
    end
    %Plot average variability (lightly)on right hand axis
    axes_sel=[1:syltypes]*2;
    subplot(syltypes,2,axes_sel);
    sim = (zero_sim-SylAveAll(:,k+2)).*100./zero_sim;
    var = (100-sim).*100./50;
    plot(SylAveAll(:,1), var, 'or');
    xlim([35 max(SylSimAll(1,1))+15]);
    ylim([0 100]);
    ylabel('Variability (ADU)');
    xlabel('Days Post Hatch');
    hold on;

    %Fit a line to the daily average variability
    fitline = polyfit(SylAveAll(:,1), var,1);
    a1 = fitline(2) % y-intercept of the fitted line
    a2 = fitline(1) % slope of fitted lines
    fit = a1+a2*SylAveAll(:,1);
    subplot(syltypes,2,axes_sel);
    plot(SylAveAll(:,1),fit,'k')
    hold off
end