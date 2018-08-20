function [PCA_set] = getPCA(DataSet)

[coeff score latent] = princomp(DataSet);
       
% for t=1:length(latent)
%     vartot=sum(latent(1:t))/length(latent);
% end

vartot = cumsum(latent)./sum(latent);

num_pc = find(vartot>.95, 1); %Determine how many prin components are required for 95% variation
%TestFeaturesPCA = score(1:num_test_syl,1:num_pc);
%TestFeaturesPCA = TestFeaturesPCA';
PCA_set = score(:,1:num_pc); %Reduce working data to only those components
end