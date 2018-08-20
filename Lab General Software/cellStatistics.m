function cellStatistics(allCells,cellArray)

%display mean of ifr as a function of correlation
for i=1:length(cellArray)
   
    yy = smooth(allCells{1,cellArray(i)}.ISINoDrug{1,1},0.1,'rlowess');
    [C,I]=max(yy(20:end));
    meanNoDrug(i)=(I+20)*4-2+rand*6;
    
    above400(i)=sum(allCells{1,cellArray(i)}.ISINoDrug{1,1}(100:200))/sum(allCells{1,cellArray(i)}.ISINoDrug{1,1}(20:200));
    
    %meanNoDrug(i)=sum(allCells{1,cellArray(i)}.ISINoDrug{1,1}.*allCells{1,cellArray(i)}.ISINoDrug{1,2});
    %activity index
    term1=sum(allCells{1,cellArray(i)}.ISINoDrug{1,1} * allCells{1,cellArray(i)}.ISINoDrug{1,2}');
    term1=term1^2;
    for j=1:length(allCells{1,cellArray(i)}.ISINoDrug{1,1})
            temp(j)=(allCells{1,cellArray(i)}.ISINoDrug{1,2}(j)^2*allCells{1,cellArray(i)}.ISINoDrug{1,1}(j));
   end
   term2=sum(temp);
          
    activityIndexNoDrug(i)=term1/term2;
    %meanDrug(i)=sum(allCells{1,cellArray(i)}.ISINoDrug{1,1}.*allCells{1,cellArray(i)}.ISINoDrug{1,2});
    corr(i)=allCells{1,cellArray(i)}.motifCorrNoDrug;
    var(i)=(100-allCells{1,cellArray(i)}.motifVarNoDrug)/50;
    age(i)=allCells{1,cellArray(i)}.dph+(allCells{1,cellArray(i)}.start-9)/8;
    %varDrug(i)=(100-allCells{1,cellArray(i)}.motifVarDrug)/50;
    %corrDrug(i)=allCells{1,cellArray(i)}.motifCorrDrug;
    %activeDrug_100(i)=allCells{1,cellArray(i)}.activeDrug{1,1}(25);
    activeNoDrug_100(i)=allCells{1,cellArray(i)}.activeNoDrug{1,1}(25);
    %activeDrug_500(i)=allCells{1,cellArray(i)}.activeDrug{1,1}(125);
    activeNoDrug_500(i)=allCells{1,cellArray(i)}.activeNoDrug{1,1}(125);
    %activeDrug_500_100(i)=allCells{1,cellArray(i)}.activeDrug{1,1}(125)/allCells{1,cellArray(i)}.activeDrug{1,1}(25);
    activeNoDrug_500_100(i)=allCells{1,cellArray(i)}.activeNoDrug{1,1}(125)/allCells{1,cellArray(i)}.activeNoDrug{1,1}(25);
%    vic(i)=allCells{1,cellArray(i)}.motifVicNoDrug;
    %vicDrug(i)=allCells{1,cellArray(i)}.motifVicDrug;
    depth(i)=allCells{1,cellArray(i)}.depth;
end
Y=activityIndexNoDrug;
%Y2=corr;
%  Y=activeDrug_500_100;
%  Y2=activeNoDrug_500_100;
% Y=varDrug;
% Y2=var;

X=age;
%Y2=varDrug;
%  fitdata = fit(X',Y','poly1')

% fitdata(1:end)
% for i=1:39
% fitY(i) = 0.007468*X(i)-0.3801;
% end
% R = corrcoef(X(1:39),Y(1:39))
  figure;
% plot(X(1:39),Y(1:39),'or',X(40:51),Y(40:51),'ok');%, X(1:39),fitY,'b');
plot(X,Y,'or');
hold on;
% plot(fitdata);
%  for i=1:length(X)
%      line([X(i) X(i)] , [Y(i) Y2(i)],'color','k');
%  end
% [b,stats] = robustfit(X,Y);
% line
