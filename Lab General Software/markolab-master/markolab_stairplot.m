function AX=sfield_stairplot(N,BINS,varargin)
%%%% generates a stair plot given histogram data
%
%
%
%


color='r';
linewidth=1;
method='l';
facecolor='r';
edgecolor='k';

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'color'
			color=varargin{i+1};
		case 'linewidth'
			linewidth=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'edgecolor'
			edgecolor=varargin{i+1};
		case 'facecolor'
			facecolor=varargin{i+1};
	end
end

% define a polygon to plot as a patch

yvec=[];
xvec=[];

for i=1:length(BINS)-1
	xvec=[xvec BINS(i:i+1)];
	yvec=[yvec repmat(N(i),[1 2])];
end

switch lower(method(1))
	case 'l'
		AX=plot(xvec,yvec,'k-','color',color,'linewidth',linewidth);
	case 'a'
		AX=area(xvec,yvec,'edgecolor',edgecolor,'facecolor',facecolor);
    case 'p'
        xvec=[-eps xvec xvec(end)+eps];
        yvec=[0 yvec 0];
        AX=patch(xvec,yvec,ones(size(xvec)),'edgecolor',edgecolor,'facecolor',facecolor);
	otherwise
		error('Did not understand plotting method');
end
