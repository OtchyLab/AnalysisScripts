function [meanDaytime, meanNighttime] = daynightStretch(cData, lData)
%Calculate the actual change in day/night stretch, accounting for the usual circadian stretch found in a control set of data.

%Set up constants for the run
sampNum = 25;
lesionTime = 0.5; %approximate time of lesion in days (i.e., hour = 24*lesionTime)
%%%%%%%%%%%%%%%%%%%%%%%%%
%Construct the stretches over the controls
%%%%%%%%%%%%%%%%%%%%%%%%%
timeSpan = 0:ceil(cData.time(end));

dn_mDur = [];
dn_mTime = [];
nd_mDur = [];
nd_mTime = [];
for i = 1:(length(timeSpan)-1)
    %%%%%%%%%
    %Day-Night Cycle
    %%%%%%%%%
    %Select out the day
    dayIdx = (cData.time > timeSpan(i)) & (cData.time < timeSpan(i+1)); 
    
    %Parse the times and lengths
    tod = cData.time(dayIdx)-timeSpan(i);
    lod = cData.int(dayIdx);
    
    if ~isempty(lod)
        dn_mDur = [dn_mDur; mean(lod(1:sampNum)), mean(lod((end-sampNum):end))];
        dn_mTime = [dn_mTime; mean(tod(1:sampNum)), mean(tod((end-sampNum):end))];
    end
    
    %%%%%%%%%
    %Day-Night Cycle
    %%%%%%%%%
        %Select out the day
    predayIdx = (cData.time > timeSpan(i)) & (cData.time < timeSpan(i+1)); 
    postdayIdx = (cData.time > (timeSpan(i)+1)) & (cData.time < (timeSpan(i+1)+1));
    
    %Parse the times and lengths
    topred = cData.time(predayIdx)-timeSpan(i);
    lopred = cData.int(predayIdx);
    topostd = cData.time(postdayIdx)-(timeSpan(i)+1);
    lopostd = cData.int(postdayIdx);
    
    if ~isempty(lopred) && ~isempty(lopostd)
        nd_mDur = [nd_mDur; mean(lopred((end-sampNum):end)), mean(lopostd(1:sampNum))];
        nd_mTime = [nd_mTime; mean(topred((end-sampNum):end)), mean(topostd(1:sampNum))];
    end

end
dn_mTime = dn_mTime.*24; %convert to hours
nd_mTime = nd_mTime.*24; %convert to hours

%Normalize to the morning length
for i = 1:size(dn_mDur,1)
    dn_mDur_norm(i,:) = dn_mDur(i,:)./dn_mDur(i,1);
end

%Limit to only those days in which rcordings span for at least 6 hours
screenIdx = (dn_mTime(:,2)-dn_mTime(:,1)) > 6;
dn_mDur_norm = dn_mDur_norm(screenIdx,:);
dn_mTime = dn_mTime(screenIdx,:);

%Normalize to the night length
for i = 1:size(nd_mDur,1)
    nd_mDur_norm(i,:) = nd_mDur(i,:)./nd_mDur(i,1);
end

%Calculate the mean daytime/nighttime adjustments
adjDaytime = mean(dn_mDur_norm(:,2),1)-mean(dn_mDur_norm(:,1),1);
adjNighttime = mean(nd_mDur_norm(:,2),1)-mean(nd_mDur_norm(:,1),1);


%%%%%%%%%%%%%%%%%%%%%%%%%
%Construct the stretches over the experimental days
%%%%%%%%%%%%%%%%%%%%%%%%%
timeSpan = 0:ceil(lData.time(end));

dnl_mDur = [];
dnl_mTime = [];
ndl_mDur = [];
ndl_mTime = [];
for i = 1:(length(timeSpan)-1)
    %%%%%%%%%
    %Day-Night Cycle
    %%%%%%%%%
    %Select out the day
    if i ~= 1
        dayIdx = (lData.time > timeSpan(i)) & (lData.time < timeSpan(i+1)); 
    else
        dayIdx = (lData.time > lesionTime) & (lData.time < timeSpan(i+1)); %On the lesion day, start after lesion
    end
    
    %Parse the times and lengths
    tod = lData.time(dayIdx)-timeSpan(i);
    lod = lData.int(dayIdx);
    
    if ~isempty(lod)
        dnl_mDur = [dnl_mDur; mean(lod(1:sampNum)), mean(lod((end-sampNum):end))];
        dnl_mTime = [dnl_mTime; mean(tod(1:sampNum)), mean(tod((end-sampNum):end))];
    end
    
    %%%%%%%%%
    %Day-Night Cycle
    %%%%%%%%%
        %Select out the day
    predayIdx = (lData.time > timeSpan(i)) & (lData.time < timeSpan(i+1)); 
    postdayIdx = (lData.time > (timeSpan(i)+1)) & (lData.time < (timeSpan(i+1)+1));
    
    %Parse the times and lengths
    topred = lData.time(predayIdx)-timeSpan(i);
    lopred = lData.int(predayIdx);
    topostd = lData.time(postdayIdx)-(timeSpan(i)+1);
    lopostd = lData.int(postdayIdx);
    
    if ~isempty(lopred) && ~isempty(lopostd)
        ndl_mDur = [ndl_mDur; mean(lopred((end-sampNum):end)), mean(lopostd(1:sampNum))];
        ndl_mTime = [ndl_mTime; mean(topred((end-sampNum):end)), mean(topostd(1:sampNum))];
    end

end
dnl_mTime = dnl_mTime.*24; %convert to hours
ndl_mTime = ndl_mTime.*24; %convert to hours

%Normalize all lengths on the basis of the mean pre-lesion length
normVal = mean(lData.int(1:sampNum));

%Normalize to the morning length
dnl_mDur_norm = dnl_mDur./normVal;

%Normalize to the night length
ndl_mDur_norm = ndl_mDur./normVal;

%Calc the deltas for the intervals
delta_dnl =dnl_mDur_norm(1:3,2) - dnl_mDur_norm(1:3,1);
delta_ndl = ndl_mDur_norm(1:2,2) - ndl_mDur_norm(1:2,1);

meanDaytime = mean(delta_dnl) - adjDaytime; %Take the mean across 3 intervals and subtract the adjustment
meanNighttime = mean(delta_ndl) - adjNighttime; %Take the mean across 2 intervals and subtract the adjustment




