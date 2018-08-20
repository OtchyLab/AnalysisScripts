function H=markolab_shadeplot(XDATA,YDATA,FACECOLOR,EDGECOLOR,METHOD)
%
%
%

if nargin<6
	METHOD=1;
else
	METHOD=2;
end

xdata=XDATA(:)';

xdata=[xdata fliplr(xdata)];
ydata=[YDATA(1,:) fliplr(YDATA(2,:))];

if METHOD==1
	H=patch(xdata,ydata,1,'facecolor',FACECOLOR,'edgecolor',EDGECOLOR);
else
	H=fill(xdata,ydata,1,'facecolor',FACECOLOR,'edgecolor',EDGECOLOR);
end
%% 
