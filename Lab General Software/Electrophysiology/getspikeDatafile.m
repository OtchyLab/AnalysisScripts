function filename = getspikeDatafile(exper, num, chan)
%This function assumes the dir command returns the filenames in alphebetical order.

filename = '';
d = dir([exper.dir,'times_',exper.birdname,'_d*chan',num2str(chan),'*']);
if(num < length(d))
    filename = d(num).name;
    currNum = extractSpikefileNumber(exper, filename);
 
else
    filename = '';
    currNum = 0;
end

if(num ~= currNum)
    lf = 1;
    rt = length(d);    
    while(true)
        mid = floor((lf + rt)/2);
        filename = d(mid).name;
        currNum = extractSpikefileNumber(exper, filename);
        if(num == currNum)
            filename = d(mid).name;
            break;
        elseif(num > currNum)
            lf = mid+1;
        elseif(num < currNum)
            rt = mid-1;
        end
        
        if(lf>rt)
            filename='0';
            return;
            %error(['getExperDatafile failed:  Filenum does not exist, or dir not returning alphebetical list.']);           
        end      
    end
end