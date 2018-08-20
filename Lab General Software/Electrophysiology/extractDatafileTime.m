function strTimeCreated = extractDatafileTime(exper,filename)

name = filename(length(exper.birdname) + length('_d') + 8:end);
strTimeCreated = name(1:15);