% projection of a vector onto a plane
function v = planeProj(plane, vec)
    v = vec - (dot(plane, vec)/(norm(plane)*norm(plane)))*plane;
end
% this isn't working either ughghghghesliuhlirueshrgelhjkgdshjkldfghkljee