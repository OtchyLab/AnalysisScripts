function vec = xRotation(coords, angle)
    % pass angle in degrees, returns new vector in the xz plane
    angle = angle * (pi/180);
    rotVec = [1, 0, 0; 0, cos(angle), -sin(angle); 0, sin(angle), cos(angle)];
    vec = rotVec * coords;
end