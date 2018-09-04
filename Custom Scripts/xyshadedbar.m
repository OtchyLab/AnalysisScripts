function [innerx, innery, outerx, outery] = xyshadedbar(x,y,dist)
%function [xs1, ys1, xs2, ys2] = xyshadedbar(x,y,dist)
% This function generates the x- and y-coordinates defining two lines that
% are equidistant from the points defined in the curve defined by (x,y).
%
% INPUTS:
%   x     --  an m x 1 array defining the x-coordinates of the central curve
%   y     --  an m x 1 array defining the y-coordinates of the central curve
%   dist  --  an m x 1 array defining the orthogonal distance from each
%             (x(i),y(i)) point. If a scalar value, the smae distance is
%             used for all points
%
% OUTPUTS:
%   xs1   -- an m x 1 array defining the x-coordinates of the outer curve offset from the central 
%   ys1   -- an m x 1 array defining the y-coordinates of the outer curve offset from the central 
%   xs2   -- an m x 1 array defining the x-coordinates of the inner curve offset from the central 
%   ys2   -- an m x 1 array defining the y-coordinates of the inner curve offset from the central 
%
% Created by TMO 8/29/2018

%Debugging switch

bDebug = false;

if bDebug
    figure(20); clf
    
    %Create a test object
    xCent = 0; yCent = 0;
    r = 1;
    
    %[x, y] = circle(xCent, yCent, r);
    % x = (rand(1,numel(x))-0.5)./10 + x;
    % y = (rand(1,numel(y))-0.5)./10 + y;
    scatter(x,y,'.'); hold on
    xlim([-1.5, 1.5]); ylim([-1.5, 1.5]);
    axis square
    
    %Test distance
    %dist = 0.25*ones(size(x));
    %dist = linspace(0, 0.25, numel(x));
end

%If dist is a scalar, use it uniformly for all points
if numel(dist) == 1
    dist = dist*ones(size(x));
end

%Estimate curve center
xc = mean(x);
yc = mean(y);

%Linear iterpolate to increase points
sc = 100;
% xHD = interp1(1:numel(x), x, linspace(1,numel(y),numel(y)*(sc+1)));
% yHD = interp1(1:numel(y), y, linspace(1,numel(y),numel(y)*(sc+1)));
xHD = interp(x, sc);
yHD = interp(y, sc);
idx = 1:sc:numel(xHD);

%Delta
diffx = diff(xHD);
diffy = diff(yHD);

%Calculate the slot at the indexed points
rg = 5;
for i = 1:numel(idx)
    start = idx(i)-rg; start = max([1,start]);
    stop = idx(i)+rg; stop = min([numel(diffy),stop]);
    spn = start:stop;
    
    dx = mean(diffx(spn));
    dy = mean(diffy(spn));
    
    
    %Simplify
    xo = xHD(idx(i)); yo = yHD(idx(i)); %original points
    x1 = xo+10e10; %test point
    
    %Solve the linear equations for the slope and y-intercept of tangent
    %and orthogonal lines
    slopes(i) = dy/dx;
    b(i) = yo - xo.*slopes(i); %tangent
    bo(i) = yo - xo.*(-1/slopes(i)); %ortho
    
    %Vectorize ortho line
    v = [xo-x1, yo-x1.*(-1/slopes(i))+bo(i)];
    len = pdist([0,0;v]);
    vn = v./len; %normallize to unit
    
    %Calculate equidistant points on either side of the curve
    xs1(i) = xo+(dist(i)*vn(1));
    xs2(i) = xo-(dist(i)*vn(1));
    
    ys1(i) = yo+(dist(i)*vn(2));
    ys2(i) = yo-(dist(i)*vn(2));
    
    %Sort into inner and outer curves
    len1 = pdist([xc,yc; xs1(i), ys1(i)]);
    len2 = pdist([xc,yc; xs2(i), ys2(i)]);
    if len1 < len2
        innerx(i) = xs1(i);
        innery(i) = ys1(i);
        outerx(i) = xs2(i);
        outery(i) = ys2(i);
    else
        innerx(i) = xs2(i);
        innery(i) = ys2(i);
        outerx(i) = xs1(i);
        outery(i) = ys1(i);
    end
    
    %Debug plotting
    if bDebug
        %line([xo, x1], [yo, x1.*slopes(i)+b(i)], 'Color', 'k')
        %line([xo, x1], [yo, x1.*(-1/slopes(i))+bo(i)], 'Color', 'g')
        line([xo, xo+(dist(i)*vn(1))], [yo, yo+(dist(i)*vn(2))], 'Color', 'y')
        line([xo, xo-(dist(i)*vn(1))], [yo, yo-(dist(i)*vn(2))], 'Color', 'y')
    end
end

%Fix occasional crossovers
for i = 1:(numel(innerx)-1)
    len1 = pdist([innerx(i), innery(i); innerx(i+1), innery(i+1)]);
    len2 = pdist([innerx(i), innery(i); outerx(i+1), outery(i+1)]);
    
    if len2<len1
        a = innerx(i+1); b = innery(i+1);
        innerx(i+1) = outerx(i+1);
        innery(i+1) = outery(i+1);
        outerx(i+1) = a;
        outery(i+1) = b;
    end
    
end

%Debug plotting
if bDebug
    scatter(xHD(idx), yHD(idx), 'ro')
    %scatter(xs1, ys1, 'ob');
    scatter(innerx, innery, 'ok');
    scatter(outerx, outery, 'ob');
    %scatter(xs2, ys2, 'ob');
end


%Test input for debugging
function [xunit, yunit] = circle(x,y,r)
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;

h = plot(xunit, yunit);
