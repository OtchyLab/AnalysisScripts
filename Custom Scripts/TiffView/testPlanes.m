% Test Planes
% create test planes to test angleAdjustment program on:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 1: -3 degree tilt from the x axis (xz planar tilt, should get a V
% angle, no U angle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get normal vector for the plane

% normal vector of a completely flat plane would be <0, 0, 1>
% we need a vector in the XZ plane that makes a 3 degree angle with said
% vector (make it a unit vector for simplicity)

nFlat = [0, 0, 1];
xyPlane = [0 0 1]; % (same as nFlat) 
xzPlane = [0 1 0];
yzPlane = [1 0 0];

nV_3 = [0.0523360, 0, 0.9986];


% Check that this gives us 3 degrees
dotProd = dot(nFlat, nV_3);
mag1 = norm(nFlat);
mag2 = norm(nV_3);
cosVAng = (dotProd)/(mag1*mag2);
Vangle = acos(cosVAng);
degVAngle = (180/3.1415926)*Vangle; % yay! 3.0001

% Now, see what angleAdjustments algorithm would return for a plane with
% normal vector nV_3 (check that there is no mistake in the latter part of
% the alg)

targVPlane = nV_3;

% V direction: (left right)
Vline = cross(targVPlane, xzPlane);
xaxis = [1 0 0];
Vangle = angleBetweenVectors(Vline, xaxis);


% U direction: (up down)
Uline = cross(targVPlane, yzPlane);
yaxis = [0 1 0];
Uangle = angleBetweenVectors(Uline, yaxis);

% degree conversions
Uangle = (360 * Uangle)/(2*pi); % 179.651 degrees (.3490 degrees)
if (Uangle > 90)
    Uangle = 180 - Uangle;
    Uangle = -Uangle;
end

Vangle = (360 * Vangle)/(2*pi); %
if (Vangle > 90)
    Vangle = 180 - Vangle; 
    Vangle = -Vangle;
end

% This seems to check out... positive x & y vals for normal vector means
% the plane would be tilted "down" into the xy plane so angle should be
% negative

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 2: 5 degree tilt from the y axis (yz planar tilt, should get a U
% angle, no V angle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nU5 = [0, .08715574, .996194698];

% Check that this gives us 5 degrees
dotProd = dot(nFlat, nU5);
Umag1 = norm(nFlat);
Umag2 = norm(nU5);
cosUAng = (dotProd)/(Umag1*Umag2);
Uangle = acos(cosUAng);
degUAngle = (180/3.1415926)*Uangle; % yay! 175.0000

% Now, see what angleAdjustments algorithm would return for a plane with
% normal vector nU5 (check that there is no mistake in the latter part of
% the alg)

targUPlane = nU5;

% V direction: (left right)
Vline = cross(targUPlane, xzPlane);
xaxis = [1 0 0];
Vangle = angleBetweenVectors(Vline, xaxis);


% U direction: (up down)
Uline = cross(targUPlane, yzPlane);
yaxis = [0 1 0];
Uangle = angleBetweenVectors(Uline, yaxis);

% degree conversions
Uangle = (360 * Uangle)/(2*pi); % 179.651 degrees (.3490 degrees)
if (Uangle > 90)
    Uangle = 180 - Uangle;
    Uangle = -Uangle;
end

Vangle = (360 * Vangle)/(2*pi); %
if (Vangle > 90)
    Vangle = 180 - Vangle; 
    Vangle = -Vangle;
end

% checks out, getting a positive angle from a normal vector in the positive
% quadrant of the yz plane... weird-ish? Probably correct based on the
% tests with real planes 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 3: 5 degree tilt from the y axis, -3 degree tilt from x axis (just
% adding the vectors from the previous two tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nUV = nV_3 + nU5;

% NORMALIZE nUV

nUV = nUV/norm(nUV);


targUVPlane = nUV;

% V direction: (left right)
Vline = planeProj(xzPlane, targUVPlane);
% PROJECTION ONTO XZ PLANE
xaxis = [1 0 0];
Vangle = angleBetweenVectors(Vline, xaxis);


% U direction: (up down)
Uline = planeProj(yzPlane, targUVPlane);
yaxis = [0 1 0];
Uangle = angleBetweenVectors(Uline, yaxis);

% degree conversions
Uangle = (360 * Uangle)/(2*pi); % 179.651 degrees (.3490 degrees)
if (Uangle > 90)
    Uangle = 180 - Uangle;
    Uangle = -Uangle;
end

Vangle = (360 * Vangle)/(2*pi); %
if (Vangle > 90)
    Vangle = 180 - Vangle; 
    Vangle = -Vangle;
end

% angles are incorrect. There must be a flaw in your method for breaking
% down the angles. Try projections?
