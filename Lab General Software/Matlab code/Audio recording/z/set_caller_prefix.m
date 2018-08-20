 function y=set_caller_prefix(caller,varh)
if strcmp(caller,'Feedback')
    if varh
        y='F'; % feedback
    else
        y='f'; % motor
    end
elseif strcmp(caller,'Playback')
    switch varh
        case 1
            y='F';  % forward
        case 0
            y='R';  % reverse
    end
end