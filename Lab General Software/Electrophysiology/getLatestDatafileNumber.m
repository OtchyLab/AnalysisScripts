function num = getLatestDatafileNum(exper)
num = 0;
d = dir([exper.dir,exper.birdname,'_d*']);
for(i = 1:length(d))
    num = max(extractDatafileNumber(exper, d(i).name),num);
end