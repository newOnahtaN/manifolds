function [ x,y,t,m ,convInd] = fbeuler_gs( dx, dy, x0, T, N, tol, maxiter)
% forward-backward euler method.
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

convInd = zeros(1,maxiter);

% computation
m = 1;
while m<maxiter && max([ norm(y-ylast), norm(x-xlast) ]) > tol

    % find x
    xlast = x;
    for i = 1:length(t)-1
        x(i+1) = x(i) + h * dx(t,x(i),ylast(i));
    end

    %find y
    ylast = y;
    for i = length(t):-1:2
        y(i-1) = y(i) - h * dy(t,xlast(i),y(i));
    end

    
    % print convergence
    %fprintf(num2str(max([ norm(y-ylast), norm(x-xlast) ])));
     %   fprintf('\n')
     convInd(m) = max([ norm(y-ylast), norm(x-xlast) ]);
     m = m+1;
end

convInd = convInd(1:m-1);
m = m-1;
    
end
