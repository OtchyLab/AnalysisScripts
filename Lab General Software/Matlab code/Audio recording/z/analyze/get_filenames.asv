%songtype_prefixA=set_songtype_prefix(caller,directedA);
%caller_prefixA=set_caller_prefix(callerA,do_feedbackA);
fnamesA=get_fnames(loadpath,birdnameA,dirnameA,'',char(callerA.songtype),char(callerA.prefix),bytes);
if length(fnamesA)>0
    namehelpsA=sort(fnamesA);
     filecounts={};
    for i=1:length(namehelpsA)
        namehelpsA{i}=namehelpsA{i}(1:end-4);
        filecounts{i}=namehelpsA{i}(end-3:end);
    end 
    namehelpA=namehelpsA{end};
    filecountA=length(namehelpsA);
else
    namehelpsA='';    namehelpA='';
    filecounts='';    filecountA=[];
end