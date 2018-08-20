function [MERGE]=frontend_collate_idxs(IMS,varargin)
%
%
%
%



low=0;
high=1;
saturation=1;
smoothing=[];
mapping=[];

nparams=length(varargin);

if mod(nparams,2)>0
	error('nndetector:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'low'
      low=varargin{i+1};
		case 'high'
			high=varargin{i+1};
		case 'saturation'
			saturation=varargin{i+1};
		case 'smoothing'
			smoothing=varargin{i+1};
		case 'mapping'
			mapping=varargin{i+1};
	end
end

if isempty(mapping)
	for i=1:length(IMS)
		mapping{i}=i;
	end
end

[rows,columns]=size(IMS{1});

for i=1:length(IMS)
	if ~all(size(IMS{i})==[rows columns])
		error('Image dimensions inconsistent');
	end
end

MERGE=zeros(rows,columns,3);

% scale images

for i=1:length(IMS)

	% scale from 0-1

	IMS{i}=(IMS{i}-min(IMS{i}(:)))./(max(IMS{i}(:))-min(IMS{i}(:)));

	% clipping, rescale

	IMS{i}=min(IMS{i},high);
	IMS{i}=max(IMS{i},low);
	IMS{i}=IMS{i}-low;
	IMS{i}=IMS{i}./(high-low);

	% map to appropriate channels

	MERGE(:,:,mapping{i})=repmat(IMS{i},[1 1 length(mapping{i})]);

end

MERGE=min(MERGE.*saturation,1);

if ~isempty(smoothing)
	MERGE=imfilter(MERGE,fspecial('disk',smoothing));
end
