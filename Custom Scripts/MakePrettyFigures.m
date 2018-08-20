function MakePrettyFigures(string)

[file, folder]=uigetfile('*.mat', 'Feature File');
vars = ['TPur360';'EPur409';'EPur433';'EPur435';'CPur407';'CPur440';'CPur432';'EPur437';'TPur357'];
for i=1:length(vars)
    load([folder file],vars(i,:));
end

figure;
eval(['plot(TPur360.' string 'x,TPur360.' string 'y,''DisplayName'',''Tutor (Pur360)'')']);
hold on;
eval(['plot(EPur409.' string 'x,EPur409.' string 'y,''r'',''DisplayName'',''Experimental (Pur409)'')']);
eval(['plot(EPur433.' string 'x,EPur433.' string 'y,''r'',''DisplayName'',''Experimental (Pur433)'')']);
eval(['plot(EPur435.' string 'x,EPur435.' string 'y,''r'',''DisplayName'',''Experimental (Pur435)'')']);
eval(['plot(CPur407.' string 'x,CPur407.' string 'y,''g'',''DisplayName'',''Control (Pur407)'')']);
eval(['plot(CPur440.' string 'x,CPur440.' string 'y,''g'',''DisplayName'',''Control (Pur440)'')']);
eval(['plot(CPur432.' string 'x,CPur432.' string 'y,''g'',''DisplayName'',''Control (Pur432)'')']);
eval(['plot(EPur437.' string 'x,EPur437.' string 'y,''r'',''DisplayName'',''Experimental (Pur437)'')']);
eval(['plot(TPur357.' string 'x,TPur357.' string 'y,''DisplayName'',''Tutor (Pur357)'')']);
end