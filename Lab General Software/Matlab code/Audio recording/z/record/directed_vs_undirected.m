songtype_prefix=set_songtype_prefix(caller,directed);
fnames=get_fnames(savepath,birdname,dirname,'',songtype_prefix,'',bytes);
if length(fnames)>0
    namehelp=char(fnames{end}(1:end-4));
    namehelp=NextName(namehelp);
    filecount=namehelp(end-3:end);
else
    filecount='0001';
end
set(hDirected,'Value',directed);
set(hUndirected,'Value',~directed);
set(hFilecount,'string',filecount);