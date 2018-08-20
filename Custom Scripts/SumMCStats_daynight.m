%create the day/night figure for motif continuity

MC = mcSet25_186;

%Format raw data
% for i = 1:3

%Move pre-lesion values
zeroIdx = find(MC(:,1) == 0, 1, 'first');

%Output arrays
time = [-1, 0:0.5:3.5];

MCout = [];
for i = 0:3
    temp = MC(zeroIdx+i, 2);
    MCout = [MCout,temp];
    
    if i == 0
        temp = MC(zeroIdx+i, 3);
        MCout = [MCout,temp];
    end
    
    temp = MC(zeroIdx+i, 4);
    MCout = [MCout,temp];
end

MC_norm = MCout./MCout(1);

day = [MC_norm(3)-MC_norm(2), MC_norm(5)-MC_norm(4), MC_norm(7)-MC_norm(6)];
night = [MC_norm(4)-MC_norm(3), MC_norm(6)-MC_norm(5)];

day_m = mean(day);
night_m = mean(night);

% MC_summary = [];
% MC_summary = [MC_summary; day_m, night_m];

figure(1);
hold on
plot([1,2], [day_m, night_m], '-s')
errorbar([1,2], mean(MC_summary,1), std(MC_summary,1)/2, '.k')

[t.dn, p.dn] = ttest(MC_summary(:,1), MC_summary(:,2));            %p = 0.0589
[t.d, p.d] = ttest(MC_summary(:,1));                                            %p = 0.102
[t.n, p.n] = ttest(MC_summary(:,2))                                             %p = 0.0529

figure(2); clf
bar([1,2], mean(MC_summary,1), 0.3, 'FaceColor', 'none'); hold on
for i = 1:4
    q(i) = plot(1:2, MC_summary(i,:), '.k');
    set(q(i), 'Marker', 'o', 'Color', [0.5, 0.5, 0.5], 'MarkerSize', 6, 'LineStyle', 'none');
end
errorbar([1,2], mean(MC_summary,1), std(MC_summary,1)/2, '.k')
xlim([0.5,2.5]); ylim([-0.1,0.4]);
ylabel('MC Change', 'FontSize', 12)
title(['MC Day/Night Recovery'], 'FontSize', 12);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2], 'XTickLabel', [{'Day'}, {'Night'}], 'YTick', [0, 0.4]);

set(gcf, 'Units', 'Inches', 'Position', [0 0 4 2.5])

