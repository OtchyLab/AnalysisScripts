function DATA=markolab_smooth(DATA,ORDER,EDGE,TYPE)
%
% OUTPUT=markolab_smooth(DATA,ORDER) 
%
% Takes a matrix of data and smooths each column 
% with a basic moving average filter (ORDER=number
% of samples in moving average)
%

if nargin<4, TYPE='b'; end
if nargin<3, EDGE='r'; end

[nsamples,ntrials]=size(DATA);

switch lower(EDGE(1))

	case 'r'

		pad_size=[2*ORDER 2*ORDER];
		front_pad=repmat(DATA(1,:),[pad_size(1) 1]);
		end_pad=repmat(DATA(nsamples,:),[pad_size(1) 1]);

	case 'z'


		pad_size=[2*ORDER 2*ORDER];
		front_pad=zeros(pad_size(1),ntrials);
		end_pad=zeros(pad_size(2),ntrials);


	case 'n'

		front_pad=[];
		end_pad=[];
		pad_size=[0 0];

	case 'w'

	otherwise


end

DATA=[front_pad;DATA;end_pad];

switch lower(TYPE(1))
	case 'b'
		b=ones(ORDER,1)/ORDER;
		a=1;
	case 'e'
		alpha=2/(1+ORDER);
		b=alpha;
		a=[1 alpha-1];
end

DATA=filter(b,a,DATA);
DATA=DATA(pad_size(1)+1:pad_size(1)+nsamples,:);

