function [LD SCORE1 SCORE2]=spikoclust_fisher_projection(DATAMAT1,DATAMAT2,varargin)
%Computes the fisher projection of two data matrices
%
%
%

%
% perhaps include plot of stability across trials (maybe threshold?)
%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

% just in case add the hot colormap at the end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fig_num'
			fig_num=varargin{i+1};
	end
end

% compute the discriminant

mu1=mean(DATAMAT1)';
mu2=mean(DATAMAT2)';

cov1=cov(DATAMAT1,1);
cov2=cov(DATAMAT2,1);

cov1=cov1+eye(size(cov1)).*1e-20;
cov2=cov2+eye(size(cov2)).*1e-20;

% pool the covariance

LD=(cov1+cov2)\(mu2-mu1);

% project both sets of points

SCORE1=DATAMAT1*LD;
SCORE2=DATAMAT2*LD;

% between class separation over within class separation

separation=((LD'*(mu1-mu2))^2)/(LD'*(cov1+cov2)*LD);
