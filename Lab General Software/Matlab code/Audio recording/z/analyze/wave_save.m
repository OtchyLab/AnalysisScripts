norm_factor=max([max(yA) abs(min(yA))]);
yA=yA/(norm_factor+0.1);
wavedirname=[dirnameA '\wf\'];
filename_path=[loadpath birdnameA '\' wavedirname num2str(filecountA) '.wav'];
wavwrite(yA,scanrateA,filename_path)
