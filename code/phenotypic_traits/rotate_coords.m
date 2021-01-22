function [xrot, yrot] = rotate_coords(x, y, x0, y0, theta)
%This function takes in matrix data x and y. Theta is the angle that the
%axes are rotated by (counterclockwise positive) for each data point.
%Returns rotated x and y values. 
xrot = cos(theta) .* (x-x0) - sin(theta) .* (y-y0) + x0; 
yrot = sin(theta) .* (x-x0) + cos(theta) .* (y-y0) + y0;
end