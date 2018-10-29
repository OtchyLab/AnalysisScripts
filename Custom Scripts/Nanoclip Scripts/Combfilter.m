%Comb filter implementation
d = fdesign.comb('notch', 'N,BW',8,0.02);
Hd = design(d,'SystemObject',true);
fvtool(Hd);

fs = 44150; fo = 510; q = 35; bw = (fo/(fs/2))/q;
[b,a] = iircomb(round(fs/fo),bw,'notch'); % Note type flag 'notch'
fvtool(b,a);