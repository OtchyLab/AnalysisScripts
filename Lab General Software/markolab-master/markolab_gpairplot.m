function [BOX_HANDLE,MED_HANDLE,WHISK_HANDLE]=markolab_gchangeplot(DATA,GRPS,varargin)
%%%% generates a stair plot given histogram data
%
%
%
%

if nargin<3, CATS=[]; end
if nargin<2, GRPS=[]; end

color='r';
linewidth=1;
outliers=0;
adjacent=0;
fill_box=1;
med_color=[0 0 1;1 0 0];
box_color=[1 0 0;0 0 1];
med_width=1;
bunching=[];
bunching_offset=4;
feature_idx=[];
box_width=.8;

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
		case 'outliers'
			outliers=varargin{i+1};
		case 'adjacent'
			adjacent=varargin{i+1};
		case 'fill_box'
			fill_box=varargin{i+1};
		case 'feature_idx'
			feature_idx=varargin{i+1};
		case 'bunching'
			bunching=varargin{i+1};
		case 'bunching_offset'
			bunching_offset=varargin{i+1};
		case 'box_width'
			box_width=varargin{i+1};
	end
end


nfeatures=length(DATA);

% position for each cell

pos_idx=ones(1,nfeatures);
if ~isempty(bunching)
	pos_idx(1:bunching:length(pos_idx))=bunching_offset;
end

pos_idx=cumsum([pos_idx]);

% now in each cell make paired plot

plot_x=[];
plot_y=[];
plot_c=[];
plot_ci=[];
nid=unique(cat(2,GRPS{:}))

%for i=1:length(DATA)
%	plot_x=[plot_x;ones(size(DATA{i}(:)))*pos_idx(i)];
%	plot_y=[plot_y;DATA{i}(:)];
%	plot_c=[plot_c;GRPS{i}(:)];
%end

for i=1:length(DATA)
	plot_x=[plot_x;pos_idx(ones(length(nid),1)*i)'+cumsum(ones(length(nid),1)*.2)];

	% get the central tendency, ci for each grp
	
	mu=[];
	ci=[];

	for j=1:length(nid)
		mu(j)=median(DATA{i}(GRPS{i}==nid(j)));
		sig=iqr(DATA{i}(GRPS{i}==j));
		ci(j,:)=[mu(j)+sig;mu(j)-sig];
	end

	plot_y=[plot_y;mu(:)];
	plot_c=[plot_c;ones(length(nid),1)*feature_idx(i)];
	plot_ci=[plot_ci;ci];
end

plot_x
plot_y
plot_ci
figure();scatter(plot_x,abs(plot_y),30,plot_c,'filled');

%size(plot_x)
%size(plot_y)
%size(plot_c)
%figure();scatter(plot_x+rand(size(plot_x))*.2,plot_y,15,plot_c);
