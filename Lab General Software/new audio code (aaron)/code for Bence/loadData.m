function [data, time, HWChannels, startSamp, timeCreated] = loadData(exper,num,chan)
%Find datafile num in the appropriate experiment folder.  Then open it.
%Num is the number of the datafile to load.
%Chan is the HWChannel number to be loaded.

filename = getExperDatafile(exper,num,chan);
if(~strcmp(filename,''))
    timeCreated = extractDatafileTime(exper,filename);
    [HWChannels, data, time, startSamp] = daq_readDatafile([exper.dir,filename]);
else
    data = [];
    time = [];
    HWChannels = [];
    startSamp = [];
    timeCreated = '';
end