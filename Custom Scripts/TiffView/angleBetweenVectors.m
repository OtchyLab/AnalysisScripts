function angle = angleBetweenVectors(v1, v2)
    % returns the acute angle between two planes (in degrees)
    dotProd = dot(v1, v2);
    mag1 = norm(v1);
    mag2 = norm(v2);
    cosAng = (dotProd)/(mag1*mag2);
    angle = acos(cosAng);
    angle = (360 * angle)/(2*pi);
end