function hh = errorbarxy2(x, y, l,u,symbol)
%ERRORBARXY Error bar plot for both x and y.
%   ERRORBARXY(X,Y,L,U, symbol) plots the graph of vector X vs. vector Y with
%   error bars specified by the vectors L and U.  L is the errorbar
%   range for each point in X and U is the errorbar range for each
%   point in Y.  Each errorbar is drawn a distance of U(i) above and
%   U(i) below the points in Y and same for X.  The vectors X,Y,L and
%   U must all be the same length. 
% 
%   For all other purposes, this program does not work. The user
%   might consider using ERRORBAR instead.
%
%   H = ERRORBARXY(...) returns a vector of line handles.
%
%   For example,
%      x = 1:10;
%      y = sin(x);
%      f = std(y)*ones(size(x));
%      e = 0.2*ones(size(y));
%      errorbarxy(x,y,e,f,'*')
%   draws symmetric error bars.

%   L. Shure 5-17-88, 10-1-91 B.A. Jones 4-5-93
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 5.19 $  $Date: 2002/06/05 20:05:14 $

%   MODIFIED BY Bilge Demirkoz for Junior lab purposes
%   Mon Nov 25 15:49:57 EST 2002

if min(size(x))==1,
  npt = length(x);
  x = x(:);
  y = y(:);
    if nargin > 2,
        if ~isstr(l),  
            l = l(:);
        end
        if nargin > 3
            if ~isstr(u)
                u = u(:);
            end
        end
    end
else
  [npt,n] = size(x);
end

if nargin == 3
  error('This function takes 5 arguments: x,y,xerr,yerr,symbol')
end

if nargin == 4
  error('This function takes 5 arguments: x,y,xerr,yerr,symbol')
end

if nargin == 2
  error('This function takes 5 arguments: x,y,xerr,yerr,symbol')
end

u = abs(u);
l = abs(l);
    
if isstr(x) | isstr(y) | isstr(u) | isstr(l)
    error('Arguments must be numeric.')
end

if ~isequal(size(x),size(y)) | ~isequal(size(x),size(l)) | ~isequal(size(x),size(u)),
  error('The sizes of X, Y, L and U must be the same.');
end

teex = (max(x(:))-min(x(:)))/100;  % make tee .02 x-distance for error bars
xl = x - teex;
xr = x + teex;
ytop = y + u;
ybot = y - u;
n = size(y,2);

% Plot graph and bars
hold_state = ishold;
cax = newplot;
next = lower(get(cax,'NextPlot'));

% build up nan-separated vector for bars
xb = zeros(npt*9,n);
xb(1:9:end,:) = x;
xb(2:9:end,:) = x;
xb(3:9:end,:) = NaN;
xb(4:9:end,:) = xl;
xb(5:9:end,:) = xr;
xb(6:9:end,:) = NaN;
xb(7:9:end,:) = xl;
xb(8:9:end,:) = xr;
xb(9:9:end,:) = NaN;

yb = zeros(npt*9,n);
yb(1:9:end,:) = ytop;
yb(2:9:end,:) = ybot;
yb(3:9:end,:) = NaN;
yb(4:9:end,:) = ytop;
yb(5:9:end,:) = ytop;
yb(6:9:end,:) = NaN;
yb(7:9:end,:) = ybot;
yb(8:9:end,:) = ybot;
yb(9:9:end,:) = NaN;

[ls,col,mark,msg] = colstyle(symbol); if ~isempty(msg), error(msg); end
symbol = [ls mark col]; % Use marker only on data part
esymbol = ['-' col]; % Make sure bars are solid

h = plot(xb,yb,esymbol); hold on

%%%%%%%%%%%%%%%%%
%% The X-errorbars %%%

teey = (max(y(:))-min(y(:)))/100;  % make tee .02 x-distance for error bars
yl = y - teey;
yr = y + teey;
xtop = x + l;
xbot = x - l;
n = size(x,2);

% Plot graph and bars
hold_state = ishold;

% build up nan-separated vector for bars
yb = zeros(npt*9,n);
yb(1:9:end,:) = y;
yb(2:9:end,:) = y;
yb(3:9:end,:) = NaN;
yb(4:9:end,:) = yl;
yb(5:9:end,:) = yr;
yb(6:9:end,:) = NaN;
yb(7:9:end,:) = yl;
yb(8:9:end,:) = yr;
yb(9:9:end,:) = NaN;

xb = zeros(npt*9,n);
xb(1:9:end,:) = xtop;
xb(2:9:end,:) = xbot;
xb(3:9:end,:) = NaN;
xb(4:9:end,:) = xtop;
xb(5:9:end,:) = xtop;
xb(6:9:end,:) = NaN;
xb(7:9:end,:) = xbot;
xb(8:9:end,:) = xbot;
xb(9:9:end,:) = NaN;

h = plot(xb,yb,esymbol); hold on

%%%%%%%%%%%%%%%%%%

h = [h;plot(x,y,symbol)]; 

if ~hold_state, hold off; end

if nargout>0, hh = h; end