function [data, time, HWChannels, startSamp, timeCreated] = readDatafile(exper,num,chan)
%Find datafile num in the appropriate experiment folder.  Then open it.

filename = getExperDatafile(exper,num,chan);
timeCreated = extractDatafileTime(exper,filename);
[HWChannels, data, time, startSamp] = daq_readDatafile([exper.dir,filename]);