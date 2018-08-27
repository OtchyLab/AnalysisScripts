function angle = getVRotation(coords)
    % returns angle in degrees
    thetaV = atan(coords(2)/coords(3)); % y/z
    angle = (360 * thetaV)/(2*pi); 
end