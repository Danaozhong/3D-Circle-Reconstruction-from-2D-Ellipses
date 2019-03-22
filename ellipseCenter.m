%
%
% Function to fit an ellipse into a set of points using a least-squares
% algorithm.
%
%
%
function [x, y] = ellipseCenter(a, b, c, d, e, f)
    num = b*b - a*c;
    x = (c*d - b*e)/num;
    y = (a*e - b*d)/num;
end