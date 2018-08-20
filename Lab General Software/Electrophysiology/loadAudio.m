function audio = loadAudio(exper,num)
%Find datafile num in the appropriate experiment folder.  Then open it.

audio = loadData(exper,num,exper.audioCh);
