function [Uangle, Vangle, finalVec] = angleAdjustments(filePath, fileName)
    % function that takes a .tif file and returns the u and v angles
    % selection range is a vector containing the minimum and maximum tiff
    % numbers you want to use to calculate the angle
    % load in tif file
    fullName = strcat(filePath, '/', fileName);
    
    % Load immages into tiff
    tiff = [];

    i = 0;
    try
        while true
            i = i + 1;
            tiff(i,:,:) = imread(fullName, i); 
            %imagesc(squeeze(tiff(i,:,:)));
            % when does this loop end?
        end
    catch ME
    end
    
    % remove edges
    edgeless_tiff = tiff(:, 11:502, 11:502);
    
    dims = size(edgeless_tiff);
    numIms = dims(1);
    height = dims(2);
    width = dims(3);
    
    % convert to binary images
    for i = 1:numIms
        for j = 1:height
            for k = 1:width
                if(edgeless_tiff(i,j,k) <= 4000)
                    edgeless_tiff(i,j,k) = 0;
                else
                    edgeless_tiff(i,j,k) = 1;
                end
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% 
    
    % determine selection range
    selectionRange = determineRange(edgeless_tiff);
    
    % select only relevant tiffs
    rel_tiffs = edgeless_tiff(selectionRange(1):selectionRange(2),:,:);
    dim = size(rel_tiffs);
    numTiffs = dim(1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find points that will be used for planar regression
    % these will be filled with all the points along the light/dark boundary
    % for each image
    allx = [];
    ally = [];
    allz = [];

    % loop through each image to get the x & y points
    % WRITE FUNCTION FOR GETTING THE X AND Y POINTS!!!
    for j = 1:numTiffs
        % second end goal: struct with x & y REGRESSION vals for each image
        curImage = squeeze(rel_tiffs(j,:,:)); % get current image
        [x, y] = findBoundary(curImage); 
        if(length(x) ~=1)
            % dont call linReg if corrcoef doesn't work... annoying
            % switch to microns now?
            x = x * 1.301; % is this the right place to do this? It is likely the problem with the final angle value
            y = y * 1.301;
            [xreg, yreg] = simpleLinReg(x,y);
            allx = [allx xreg];
            ally = [ally yreg];
            zadd = zeros(1, length(xreg)) + j;
            allz = [allz zadd];
        end

        [r, m, b] = regression(x, y);
        regressStats(1,j) = r;
        regressStats(2,j) = m;
        regressStats(3,j) = b;
    end

    allz = allz(1:length(allx));
    
    % Perform planar regression using x, y, and z points 
    targPlane = planarRegression(allx,ally, allz);

    %%%%%%%%%%%%%%%%%%%%%%
    
    Vangle = getVRotation(targPlane);
    
    rotatedTarg = xRotation(targPlane, Vangle); 
    
    zAxis = [0 0 1];
    Uangle = angleBetweenVectors(rotatedTarg, zAxis);
    
    finalVec = yRotation(rotatedTarg, Uangle);
    
end
%%%%%%%%%%%%%%%%%%%%%
%   Subfunctions    %
%%%%%%%%%%%%%%%%%%%%%
function selectionRange = determineRange(tiff) 
    % go through a tiff stack and determine which images to use in angle
    % calculations based on the amount of white and dark pixels
    lowThresh = 85781;
    highThresh = 237697;
    dim = size(tiff);
    for i = 1:dim(1)
        if(sum(sum(squeeze(tiff(i,:,:)))) > lowThresh)
            selectionRange(1) = i;
            start = i+1;
            break;
        end
    end
    for j = start:dim(1)
        if(sum(sum(squeeze(tiff(j,:,:)))) > highThresh)
            selectionRange(2) = j - 1;
            break;
        end
    end
end

function [xreg, yreg] = simpleLinReg(x,y)
    xbar = mean(x);
    xstd = std(x);
    ybar = mean(y);
    ystd = std(y);
    rmat = corrcoef(x,y); % I hate matlab
    r = rmat(2,1);

    b = r * (ystd/xstd);
    a = ybar - b*xbar;

    % Y = bX + a
    % 
    xreg = x;
    yreg = b*xreg + a;
end

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

function angle = angleBetweenVectors(v1, v2)
    % returns the acute angle between two planes (in radians)
    dotProd = dot(v1, v2);
    mag1 = norm(v1);
    mag2 = norm(v2);
    cosAng = (dotProd)/(mag1*mag2);
    angle = acos(cosAng);
end

%function binIm = convertToBinary(tiff)
%
%end

