function [stimIdx, timeIdx] = createStimIdx(triggers)
%Function reads in a cell array of stimulation trigger (TTL) pulses
%(~0-5VDC). Using threshold crossings, it outputs both stim condition
%(i.e., stim/no-stim) and time index (trigger sample point).
%
%INPUT:
%   triggers = 1xN cell array (where N is the number of motifs) in which
%               each vector is a time series containing TTL signals corresponding
%               to stimulator GO commands
%OUTPUT:
%   stimIdx = 1xN boolean array where 0 = no stiml 1 = stim
%   timeIdx = 1xN cell array identifying the sample number corresponding
%               to the rising edge of a detected TTL pulse. Note that there may 
%               be more than one of these per motif (hence the use of a cell array)
%
%Written by TMO 12/12/17

%Input condition check
if isempty(triggers)
    return
end

%Predefine the output variables
stimIdx = false(size(triggers));
timeIdx = cell(size(triggers));

%Set threshold for TTL detection
thresh = 2.5; %5VDC signal, so something lower should work fine, I think...

%Cycle through the triggers and fill in the arrays
rends = numel(triggers);
for i = 1:rends
    %Shorten
    vx = triggers{i};
    
    %Threshold crossings (only look at rising edge)
    p = 2:numel(vx);
    ups = find(vx(p)>thresh & vx(p-1)<thresh);
    
    %Update indices as appropriate
    if ~isempty(ups)
        stimIdx(i) = true;
        timeIdx{i} = ups;
    end
end
    
    
    