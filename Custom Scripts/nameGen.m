function [names, titles] = nameGen(pntr, chan)

chan = num2str(chan);
if pntr == 1
    freqs = '180-500Hz';
elseif pntr == 2
     freqs = '180-1000Hz'; 
elseif pntr == 3 
    freqs = '300-1000Hz';
elseif pntr == 4
    freqs = '300-6000Hz';
elseif pntr == 5
    freqs = '500-1000Hz';
elseif pntr == 6
    freqs = '1000-6000Hz';
end

 names = {['Grn141_141210_dataset ch' chan ' ' freqs '.mat'];...
                ['Grn141_141211_dataset ch' chan ' ' freqs '.mat'];...
                ['Grn141_141214_dataset ch' chan ' ' freqs '.mat']};
                
titles = ['Grn141 HVC Activity (Ch' chan ', ' freqs ')'];