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