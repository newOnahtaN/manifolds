function [ x,y,t,m ] = fbtrap_jb( dx, dy, x0, T, N, tol, maxiter)
% forward-backward trapezoid method with a Jacobi update scheme.
%
% inputs:
%
% dx -- function of t, x, and y
% dy -- function of t, x, and y
% x0 -- initial x
% T -- boundary position
% N -- number of grid points
% tol -- tolerance for convergence
% maxiter -- Maximum number of iterations
%
% outputs: 
%
% x -- array of x values
% y -- array of y values
% t -- array of t values
% m -- number of steps to convergence
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make a grid
h = (T+1)/N;
t = 0:h:T;
x = x0*ones(1, length(t));
% initialize y as zeros. We can discuss if this is valid or not later.
y = zeros(1, length(t));
xlast = x*0;
ylast = y*0;

% computation
m = 1;
while m<maxiter && max([ norm(y-ylast), norm(x-xlast) ]) > tol

    % find
    xlast = x;
    for i = 1:length(t)-1
        x(i+1) = xlast(i) + h * dx(t(i),(xlast(i)+xlast(i+1))/2,(ylast(i)+ylast(i+1))/2);
    end

    %find y
    ylast = y;
    for i = length(t):-1:2
        y(i-1) = ylast(i) - h * dy(t(i),(xlast(i)+xlast(i-1))/2,(ylast(i)+ylast(i-1))/2);
    end

    m = m+1;
    % print convergence
    fprintf(num2str(max([ norm(y-ylast), norm(x-xlast) ])));
        fprintf('\n')
end

    
end
