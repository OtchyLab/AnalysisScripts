function [x y] = findBoundary(tiff)
    approach = approachDir(tiff);
    
    [rows, cols] = size(tiff);
    x = [];
    y = []; 
    idx = 1; 
    switch approach
        case 1 % top
            counter = 1;
            for i = 1:cols
                while ( tiff(counter, i) == 1 && counter < rows)
                    counter = counter + 1; 
                end
                if(counter ~= rows) % black pixel found
                    y(idx) = counter;
                    x(idx) = i; 
                    idx = idx + 1;
                end
                counter = 1;
            end
        case 2 % left
            counter = 1;
            for i = 1:rows
                while (tiff(i, counter) == 1 && counter < cols)
                    counter = counter + 1;
                end
                if(counter ~= cols) % black pixel found
                    y(idx) = i;
                    x(idx) = counter; 
                    idx = idx + 1;
                end
                counter = 1;
            end
        case 3 % bottom
            counter = rows;
            for i = 1:cols
                while (counter > 0 && tiff(counter,i) == 1)
                    counter = counter - 1; 
                end
                if(counter ~= 0) % black pixel found
                    y(idx) = counter;
                    x(idx) = i; % columns are x vals
                    idx = idx + 1; 
                end
                counter = rows;
            end
        otherwise % right
            counter = cols; 
            for i = 1:rows
                while (counter > 0 && tiff(i, counter) == 1)
                    counter = counter - 1; 
                end
                if(counter ~= 0) % black pixel found
                    y(idx) = i;
                    x(idx) = counter;
                    idx = idx + 1; 
                end
                counter = cols; 
            end
    end
    
    x = fliplr(x);
    
    % returns x and y :)
    
end
%%%% Subfunctions %%%
function approach = approachDir(tiffIm)
    % this function needs to give a consistently good way to find the
    % boundary between the light and dark pixels in a binary image
    
    % argument is a 2d array of 0's and 1's
    
    % evaluate the edges: which edges are completely dark?
    % the way you look at the image will depend on this 
    
    % return 1 for top, 2 for left, 3 for bottom, 4 for right
    
    [r, c] = size(tiffIm); 
    
    topWhite = 0;
    for j = 1:c
        if (tiffIm(1, j) == 1) % if the pixel is white
            topWhite = topWhite + 1; 
        end
    end
    
    mostWhite = topWhite; 
    approach = 1; % approach from top by default
    
    bottomWhite = 0;
    for j = 1:c
        if (tiffIm(r, j) == 1) % if the pixel is white
            bottomWhite = bottomWhite + 1; 
        end
    end
    
    if bottomWhite > mostWhite
        mostWhite = bottomWhite; 
        approach = 3; 
    end
    
    rightWhite = 0;
    for j = 1:r
        if (tiffIm(j, c) == 1) % if the pixel is white
            rightWhite = rightWhite + 1; 
        end
    end
    
    if rightWhite > mostWhite
        mostWhite = rightWhite; 
        approach = 4;
    end
    
    leftWhite = 0;
    for j = 1:r
        if (tiffIm(j, 1) == 1) % if the pixel is white
            leftWhite = leftWhite + 1; 
        end
    end
    
    if leftWhite > mostWhite
        mostWhite = leftWhite; 
        approach = 2;
    end
    
end