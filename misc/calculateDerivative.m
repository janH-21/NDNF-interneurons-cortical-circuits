function [dy_dx, d2y_dx2] = calculateDerivative(x,y)
% INPUTS:
%   x - 1D vector of x values of function
%   y - 1D vector of y values of funtion
% OUTPUTS:
%   dy_dx - first derivative of y with respect to x
%   d2y_dx2 - second derivative of y with respect to x
% NOTE:
% Vectors x and y need to be of equal lenght.

    % preallocate empty vectors
    dy_dx = zeros(1, length(y));
    d2y_dx2 = zeros(1, length(y));
    
    % calculate ds
    dx = gradient(x);
    dy = gradient(y);
    
    % calculate dy_dx (using a loop instead of ./ to not run into trouble with memory
    % when using very large vecotrs)
    for iter = 1:length(dx)
        dy_dx(iter) = dy(iter)/dx(iter);
    end
    
    % calculate d2y_dx2
    d2y = gradient(dy_dx);
    for iter = 1:length(dx)
        d2y_dx2(iter) = d2y(iter)/dx(iter);
    end
    
    
end