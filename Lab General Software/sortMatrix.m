function [value,index]=sortMatrix(data)
    
%sorts a correlation matrix where the minimum value is -1
for i=1:(size(data,1)*(size(data,1)-1))/2
      maxdata=max(data(:));
      value(i)=maxdata;
      [r,c]=find(data==maxdata,1,'first');
      data(r,c)=-1;
      index(i,:)=[r,c];
end
