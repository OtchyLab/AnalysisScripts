function [FeatScaled] = ScaleFeats(FeatureSet, type)
%Scales all of the paramaters in the specified matrix to the empirically
%determined median absolute difference from the mean.

%Get size of FeatureSet for future calculations
[m,n] = size(FeatureSet);

%Safety to prevent inadvertant overwrite
if nargin == 1
    type = 1;
end

if type == 666
    %Determine the mean for each paramter in in the set
    FeatMeans = mean(FeatureSet);

    %Determine the difference of each point from the mean
    FeatMeans = ones(m,1)*FeatMeans;
    FeatDeviations = abs(FeatureSet-FeatMeans);

    %Get the median of the deviations for each parameter
    FeatMedians = median(FeatDeviations);
else
    %These value determined experimentally 01/06/10
    FeatMedians = [0.0454622631715250,724.984533730285,376.720986859473,405.795758917253,319.747858391467,0.0114011114061349,31.6657937501424,10.3833824444278,5.88012484818616,25.5863709361087,33.2300569344639,31.1925666657626,29.0274723425039,1.16759727979522,1.20808605838609,1.85292219878446,0.771091956520207,0.906378072714191,1.17108614800320,1.31524525635668,0.892283336590318,44.7108183684006,37.3180199436283,40.3935302121631,30.2212452254079,49.7417406581550,37.4546097366776,28.7860188761821,74.0320422982997,8.59707336703535,11.2607469892210,14.4039849880523,10.4639504257984];
end

%Scale the original FeatureSet to the MADs from the mean
FeatScaled = FeatureSet./(ones(m,1)*FeatMedians);

end