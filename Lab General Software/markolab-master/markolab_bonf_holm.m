function [P_COR,H]=markolab_bonf_holm(P,ALPHA)
%computes the bonf holm MC correction
%%
%

if ~isvector(P)
	error('P must be a vector...');
else
	P=P(:)';
end

H=zeros(size(P));
P_COR=nan(size(P));

% sort the p-values and multiply by the correction

[sort_p,sort_idx]=sort(P,'ascend');

m=length(sort_p);
correction=m:-1:1;
tmp=sort_p.*correction;
P_COR(sort_idx)=tmp;
H=P_COR<ALPHA;


