function time = extractExperFilenumTime(exper,filenum)

ch = exper.audioCh;
filename = getExperDatafile(exper,filenum, ch);
filetimestr = extractDatafileTime(exper,filename);

y = str2num(filetimestr(1:4));
m = str2num(filetimestr(5:6));
d = str2num(filetimestr(7:8));
th = str2num(filetimestr(10:11));
tm = str2num(filetimestr(12:13));
ts = str2num(filetimestr(14:15));

time = datenum(y,m,d,th,tm,ts);
