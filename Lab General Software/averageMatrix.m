function matrixVector=averageMatrix(symmMatrix,flag)
%takes a pseudo-symmetric matrix (like those from SAP) and returns a vector where the 
%elements are ((1,2)+(2,1)/2 ((1,3)+(3,1))/2...(1,n)(2,3)...(2,n) etc. 
matrixVector=[];
for i=1:size(symmMatrix,1)-1
    for j = i+1:size(symmMatrix,1)-1
        if (flag==1)
            matrixVector(end+1)=(symmMatrix(i,j)+symmMatrix(j,1))/2;
        else
            matrixVector(end+1)=symmMatrix(i,j);
        end
    end
end
