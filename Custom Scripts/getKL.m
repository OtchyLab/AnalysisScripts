function KL = getKL(P,Q)

%Smooth using SGTS
Psm = sgts(P);
Qsm = sgts(Q);

%Determine the bounds of the subsequent analysis based on the 
t = P+Q;
x = sum(t,2); y = sum(t,1); 
minx = find(x ~=0,1,'first'); maxx = find(x ~=0,1,'last');
miny = find(y ~=0,1,'first'); maxy = find(y ~=0,1,'last');

%Reduce matrices to their minimal forms and normalize
Ps = Psm(minx:maxx, miny:maxy); Ps = Ps/sum(Ps(:));
Qs = Qsm(minx:maxx, miny:maxy); Qs = Qs/sum(Qs(:));

%Calculate KL divergence on the remaining data (in bits)
KL = sum(Ps(:).*log(Ps(:)./Qs(:)))/log(2);
