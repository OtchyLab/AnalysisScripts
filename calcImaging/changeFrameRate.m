function audioData = changeFrameRate(raw, rate_orig, rate_new)
    % Takes raw audio data at 48000
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fixing the frequency of audio recordings! :)
    %vq = interp1(x,v,xq) returns interpolated values of a 1-D function at 
    %specific query points using linear interpolation. Vector x contains the 
    %sample points, and v contains the corresponding values, v(x). Vector xq 
    %contains the coordinates of the query points.
    % assumes raw is a column vector

    x = 1:numel(raw);
    v = raw';

    fs = rate_new; %desired frequency
    n = round((numel(raw)/rate_orig) * fs);

    xq = linspace(1, numel(v), n);

    audioData = interp1(x, v, xq);

    % yayyy works!!!

end