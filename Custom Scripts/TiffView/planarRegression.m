function ans = planarRegression(xs,ys,zs)
% calculates the plane of best fit for a group of points
    xsum = sum(xs);
    ysum = sum(ys);
    zsum = sum(zs); 
    onesum = length(xs);
    xymult = sum(xs.*ys);
    xzmult = sum(xs.*zs);
    yzmult = sum(ys.*zs);
    xxmult = sum(xs.*xs);
    yymult = sum(ys.*ys);
    %zzmult = sum(zs.*zs);
    
    A = [xxmult xymult xsum; xymult, yymult, ysum; xsum ysum onesum];
    A = inv(A); 
    
    colVec = [xzmult; yzmult; zsum];
    
    ans = A*colVec;
    
end