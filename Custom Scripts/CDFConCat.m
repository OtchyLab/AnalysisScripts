function [output]=CDFConCat

[file, folder]=uigetfile('*.mat', 'Feature File');
vars = ['TPur360';'EPur409';'EPur433';'EPur435';'CPur407';'CPur440';'CPur432';'EPur437';'TPur357'];
output = [];
TempVecAM = [];
TempVecFM = [];
TempVecEnt = [];
TempVecGP = [];

for i=1:length(vars)
    load([folder file],vars(i,:));
    %eval(['TempVec = ' vars(i,:) '.ecdf.yproj.AM, ' vars(i,:) '.ecdf.yproj.Ent, ' vars(i,:) '.ecdf.yproj.GP;']);
    eval(['TempVecAM = [TempVecAM; ' vars(i,:) '.ecdf.yproj.AM''];'])
    eval(['TempVecFM = [TempVecFM; ' vars(i,:) '.ecdf.yproj.FM''];'])
    eval(['TempVecEnt = [TempVecEnt; ' vars(i,:) '.ecdf.yproj.Ent''];'])
    eval(['TempVecGP = [TempVecGP; ' vars(i,:) '.ecdf.yproj.GP''];'])
    %output = [output; TempVec'];
end

%Z-score each feature array
m_AM = mean(TempVecAM);
s_AM = std(TempVecAM);
z_AM = (TempVecAM-m_AM)./s_AM;

m_FM = mean(TempVecFM);
s_FM = std(TempVecFM);
z_FM = (TempVecFM-m_FM)./s_FM;

m_Ent = mean(TempVecEnt);
s_Ent = std(TempVecEnt);
z_Ent = (TempVecEnt-m_Ent)./s_Ent;

m_GP = mean(TempVecGP);
s_GP = std(TempVecGP);
z_GP = (TempVecGP-m_GP)./s_GP;

%Breakout to output
bo = 1:101:(length(vars)+1)*101;

output = [[z_AM(bo(1):bo(2)-1)'],[z_FM(bo(1):bo(2)-1)'],[z_Ent(bo(1):bo(2)-1)'],[z_GP(bo(1):bo(2)-1)']];
for i=2:length(vars)
    output = [output; [z_AM(bo(i):bo(i+1)-1)'],[z_FM(bo(i):bo(i+1)-1)'],[z_Ent(bo(i):bo(i+1)-1)'],[z_GP(bo(i):bo(i+1)-1)']];
end



end