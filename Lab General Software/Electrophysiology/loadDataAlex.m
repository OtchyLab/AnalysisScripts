function data_out=loadDataAlex(rec,channel)

if (rec<10)
    fill='000';
elseif (rec<100)
    fill='00';
else
    fill='0';
end
filename=['adult_0714D-f' fill num2str(rec) '.daq']; 

data = daqread(filename);
data_out=data(:,channel);
end
