function hOut = makewedge(rho1, rho2, theta1, theta2, color)
%MAKEWEDGE Plot a wedge.
% MAKEWEDGE(rho1, rho2, theta1, theta2, color) plots a polar
% wedge bounded by the given inputs. The angles are in degrees.
%
% h = MAKEWEDGE(...) returns the patch handle.

% ang = linspace(theta1/180*pi, theta2/180*pi);
ang = linspace(theta1, theta2);
[arc1.x, arc1.y] = pol2cart(ang, rho1);
[arc2.x, arc2.y] = pol2cart(ang, rho2);
x = [arc1.x arc2.x(end:-1:1)];
y = [arc1.y arc2.y(end:-1:1)];

h = patch(x, y, color);
if ~ishold
    axis equal tight;
end
if nargout > 0
     hOut = h;
end