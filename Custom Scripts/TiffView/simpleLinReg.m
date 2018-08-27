function [xreg, yreg] = simpleLinReg(x,y)
    xbar = mean(x);
    xstd = std(x);
    ybar = mean(y);
    ystd = std(y);
    rmat = corrcoef(x,y); % I hate matlab
    r = rmat(2,1);

    b = r * (ystd/xstd);
    a = ybar - b*xbar;

    % Y = bX + a
    % 
    xreg = x;
    yreg = b*xreg + a;
end