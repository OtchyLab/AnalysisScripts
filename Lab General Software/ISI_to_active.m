
function active=ISI_to_active(IFR, IFR_bins)

tim(1:201)=1./IFR_bins';
isiTimeWeighted(1,1:201)=IFR.*tim;

for k=1:201
    active(k)=sum(isiTimeWeighted(k:201))/sum(isiTimeWeighted(1:201));

end
figure;
semilogy(IFR_bins,active);