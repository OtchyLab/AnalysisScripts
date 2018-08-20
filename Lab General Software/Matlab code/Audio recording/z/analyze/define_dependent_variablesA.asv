define_callers;

bytes=8000; 
path(path,'..'); %sets path

if ~exist('callerA') 
    callerA={};
    caller_iA=get_index(callers.name,callerA_default);%1 if 'Feedback; 2 if 'Playback'
    callerA.name=callerA_default; %default now is feedback
    callerA.prefix=callers.prefix{caller_iA}(1); %(now 'f') for motor 
    callerA.prefix_text=callers.prefix_text{caller_iA}(1); 
    callerA.songtype=callers.songtype{caller_iA}(1); %now directed
    callerA.songtype_text=callers.songtype_text{caller_iA}(1);
end
if ~exist('do_plotA')
    do_plotA=do_plot_defaultA; %=1
end
if ~exist('online')
    online=online_default; %=0
end
%if ~exist('caller_prefixA')
%    caller_prefixA=callers.prefix(caller_iA;
%end
%if ~exist('songtypeA')
%    songtypeA=1;
%end
%songtype_prefixA=set_songtype_prefix(caller,directedA);
if exist('savepath')
    loadpath=savepath;
else
    loadpath='/Matlab-6.5/Data/';
end
