function y=set_songtype_prefix(caller,varh)
if strcmp(caller,'Feedback')
    if varh
        y='D';   % directed
    else
        y='U';   % undirected
    end
elseif strcmp(caller,'Playback')
    if varh
        y='B';  % BOS
    end
end