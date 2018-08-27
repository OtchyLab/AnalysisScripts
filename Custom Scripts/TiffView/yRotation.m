function vec = yRotation(coords, angle)
    % take angle in degrees
    angle = -angle * (pi/180);
    rotVec = [cos(angle) 0 sin(angle); 0 1 0; -sin(angle) 0 cos(angle)];
    vec = rotVec * coords;
end