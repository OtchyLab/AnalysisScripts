hold on; h = plotCorrelationDistribution(RRotvalues20{1}, .05); set(h,'Color','red'); set(h,'LineStyle','-'); set(h,'LineWidth',1);
hold on; h = plotCorrelationDistribution(Rvalues20{1}, .05); set(h,'Color','red'); set(h,'Marker','o'); set(h,'LineWidth',1);
hold on; h = plotCorrelationDistribution(allRRotValues20, .05); set(h,'Color','black'); set(h,'LineStyle','-'); set(h,'LineWidth',1);
hold on; h = plotCorrelationDistribution(allRValues20, .05); set(h,'Color','black'); set(h,'Marker','o'); set(h,'LineWidth',1);

 std(allRRotValues20)
 std(allRValues20)
 std(RRotvalues20{1})
 std(Rvalues20{1})

 mean(allRRotValues20)
 mean(allRValues20)
 mean(RRotvalues20{1})
 mean(Rvalues20{1})
 