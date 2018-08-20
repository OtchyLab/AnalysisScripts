function [FIGNUM,LABELS,SELECTION,DIMS]=markolab_clust_cut(FEATURES,FEATURE_LABELS,CLUSTERFUN,DISPLAY_LIM)
% general purpose GUI for cluster cutting, pass in FEATURES, returns LABELS
% TODO: subsample display

disp('Complete selection by closting the clustering window...')

if nargin<4 | isempty(DISPLAY_LIM)
	DISPLAY_LIM=inf;
end

if nargin<3 | isempty(CLUSTERFUN)
	CLUSTERFUN=@(x,n) kmeans(x,n);
end

[nsamples,nfeatures]=size(FEATURES);

LABELS=ones(nsamples,1);
SELECTION=1;
DIMS=[];

if nargin<2 | isempty(FEATURE_LABELS)
	for i=1:nfeatures
		FEATURE_LABELS{i}=[ 'Feature ' num2str(i) ];
	end
end

if nsamples>DISPLAY_LIM
	disp_idx=unique(round(linspace(1,nsamples,DISPLAY_LIM)));
	disp([sprintf('Displaying %g out of %g total points)',length(disp_idx),nsamples)]);
else
	disp_idx=1:nsamples;
end

ndims=3;

main_window=figure('Visible','off','Position',[360,500,700,600],'Name','Markolab Clust Cut','NumberTitle','off');
plot_axis=axes('Units','pixels','Position',[50,50,400,400]);

pop_up_x= uicontrol('Style','popupmenu',...
	'String',FEATURE_LABELS,...
	'Position',[400,90,90,25],'call',@change_plot,...
	'Value',1);
pop_up_x_text= uicontrol('Style','text',...
	'String','X',...
	'Position',[405,130,50,45]);


pop_up_y= uicontrol('Style','popupmenu',...
	'String',FEATURE_LABELS,...
	'Position',[495,90,90,25],'call',@change_plot,...
	'Value',min(nfeatures,2));
pop_up_y_text= uicontrol('Style','text',...
	'String','Y',...
	'Position',[500,130,50,45]);

pop_up_z= uicontrol('Style','popupmenu',...
	'String',FEATURE_LABELS,...
	'Position',[595,90,90,25],'call',@change_plot,...
	'Value',min(nfeatures,3));
pop_up_z_text= uicontrol('Style','text',...
	'String','Z',...
	'Position',[600,130,50,45]);

pop_up_clusters= uicontrol('Style','popupmenu',...
	'String',[1:9],...
	'Position',[475,210,75,25],'value',1);
pop_up_clusters_text= uicontrol('Style','text',...
	'String','N(clusters) auto-clust',...
	'Position',[500,250,100,45]);

pop_up_choice= uicontrol('Style','popupmenu',...
	'String',[1:20],...
	'Position',[475,330,75,25],'value',1,'call',@change_selection);
pop_up_choice_text= uicontrol('Style','text',...
	'String','Cluster selection',...
	'Position',[500,370,100,45]);

push_draw_mode= uicontrol('Style','pushbutton',...
	'String','Draw cluster (X and Y only)',...
	'Position',[500,475,180,35],'value',0,...
	'Call',@change_cluster);
push_auto_mode = uicontrol('Style','pushbutton',...
	'String','Auto-clust (all checked feat.)',...
	'Position',[500,425,180,35],'value',0,...
	'Call',@change_cluster);

rows=ceil(length(FEATURE_LABELS)/5);


status_text= uicontrol('Style','text',...
	'String','Status: Waiting for input...',...
	'Position',[600,550,150,45]);


nrows=ceil(length(FEATURE_LABELS)/5);
ncolumns=5;

feature_panel= uipanel('Parent',main_window,...
	'Title','Auto-clust features',...
	'Units','Pixels',...
	'Position',[50 540-nrows*35 80*ncolumns 50*nrows]);

i=1;
while i<=length(FEATURE_LABELS)
	row=ceil(i/5);
	column=mod(i,5);
	if column==0, column=5; end
	cluster_dims{i}=uicontrol('parent',feature_panel,...
		'Style','checkbox',...
		'String',FEATURE_LABELS{i},...
		'Units','Normalized',...
		'Value',0,'Position',[(column-1)*(.98/ncolumns),1-row*(1/nrows),1/ncolumns,1/nrows]);
	set(cluster_dims{i},'Units','Normalized')
	i=i+1;
end

set(cluster_dims{1},'Value',1);

% now align everything and send the main_window handle to the output

align([pop_up_clusters,pop_up_clusters_text,...
	pop_up_choice,pop_up_choice_text,push_draw_mode,push_auto_mode,status_text],'Center','None');
align([pop_up_x,pop_up_x_text],'Center','None');
align([pop_up_y,pop_up_y_text],'Center','None');
align([pop_up_z,pop_up_z_text],'Center','None');

FIGNUM=gcf;

% run change_plot, which updates the plot according to the defaults

change_plot();

set(pop_up_choice,'string',[1:length(unique(LABELS))]);
set(pop_up_choice,'value',1);

set([main_window,plot_axis,pop_up_x,pop_up_x_text,pop_up_y,pop_up_y_text,pop_up_z,...
	pop_up_z_text,pop_up_clusters,pop_up_clusters_text,pop_up_choice,pop_up_choice_text,...
	push_draw_mode,push_auto_mode,status_text,feature_panel],'Units','Normalized');

movegui(main_window,'center')
set(main_window,'Visible','On')
waitfor(main_window);

%% Callbacks

% this callback changes the plot and returns the sum of the distances
% from the centroid for each point in a cluster

function change_plot(varargin)

viewdim(1)=get(pop_up_x,'value');
viewdim(2)=get(pop_up_y,'value');
viewdim(3)=get(pop_up_z,'value');

clusters=unique(LABELS);
clusternum=length(clusters);
% clear the plot axis

cla;
colors=colormap(['lines(' num2str(clusternum) ')']);

switch ndims
	case 2
		for i=1:clusternum
			clusterid{i}=num2str(i);
			points=find(LABELS(disp_idx)==clusters(i));
			h(:,i)=plot(FEATURES(disp_idx(points),viewdim(1)),FEATURES(disp_idx(points),viewdim(2)),...
				'o','markerfacecolor',colors(i,:),'markeredgecolor','none');hold on
		end

	case 3
		for i=1:clusternum
			clusterid{i}=num2str(i);
			points=find(LABELS(disp_idx)==clusters(i));
			h(:,i)=plot3(FEATURES(disp_idx(points),viewdim(1)),FEATURES(disp_idx(points),viewdim(2)),...
			FEATURES(disp_idx(points),viewdim(3)),'o','markerfacecolor',colors(i,:),'markeredgecolor','none');hold on
		end

end

% label everything
grid on
view(ndims)

xlabel(FEATURE_LABELS{viewdim(1)});
ylabel(FEATURE_LABELS{viewdim(2)});
zlabel(FEATURE_LABELS{viewdim(3)});

L=legend(h,clusterid,'Location','NorthEastOutside');legend boxoff
set(L,'FontSize',20,'FontName','Helvetica')

end

function change_selection(varargin)

tmp=get(pop_up_choice,'string');
SELECTION=str2num(tmp(get(pop_up_choice,'value')));

end

function change_cluster(varargin)
%
%

disp('Updating clusters...');


choices=get(pop_up_clusters,'string');
clusternum=str2num(choices(get(pop_up_clusters,'value')));

viewdim(1)=get(pop_up_x,'value');
viewdim(2)=get(pop_up_y,'value');
viewdim(3)=get(pop_up_z,'value');

draw_mode=get(push_draw_mode,'value');

if draw_mode
	set(status_text,'string','Status:  Updating clusters (draw mode)...');
else
	set(status_text,'string','Status:  Updating clusters (auto-clust)...')
end
drawnow;

dim=[];
for i=1:length(cluster_dims)
	val=get(cluster_dims{i},'value');
	if val==1
		dim=[dim i];
	end
end

if isempty(dim)
	dim=viewdim;
	for i=1:length(dim)
		set(cluster_dims{dim(i)},'value',1);
	end
end

DIMS=dim;

choices=get(pop_up_clusters,'string');
clusternum=str2num(choices(get(pop_up_clusters,'value')));

% cluster, user-specified function is passed

LABELS=CLUSTERFUN(FEATURES(:,dim),clusternum);

% clear the plot axis

cla;

% plot in either 2 or 3 dims
%

if draw_mode

	drawnow;
	ndims=2;
	plot(FEATURES(disp_idx,viewdim(1)),FEATURES(disp_idx,viewdim(2)),'o','markerfacecolor','b');view(2);
	hold on
	disp('Select the corners of the enclosing polygon then press RETURN to continue...');
	hold off;

	LABELS=ones(nsamples,1);
	counter=2;

	response2=[];

	while isempty(response2)

		[xv,yv]=ginput;
		k=convhull(xv,yv);
		hold on;
		plot(xv(k),yv(k),'b-','linewidth',1.25);
		choice=inpolygon(FEATURES(:,viewdim(1)),FEATURES(:,viewdim(2)),xv(k),yv(k));
		LABELS(choice==1)=counter;
		response2=questdlg('Done drawing or keep going?','Keep going?',...
			'Done','Continue','Done');

		switch lower(response2(1))
			case 'd'
				break;
			case 'c'
				response2=[];
			otherwise
				response2=[];
		end

		counter=counter+1;


	end

	clusternum=counter;
else
	ndims=3;
	LABELS=CLUSTERFUN(FEATURES(:,dim),clusternum);
end

set(pop_up_choice,'value',1);
set(pop_up_choice,'string',[1:length(unique(LABELS))]);
change_plot();
set(status_text,'string','Status:  Clustering complete, waiting for input...');

end

end
