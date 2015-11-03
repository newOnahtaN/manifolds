clear; clc; close all
% test stability on a linear test system

% evolution rules
% define equations
dx = @(t,x,y) -x+2*y+(x-y).^2;
dy = @(t,x,y) (x-y).^2+y;
% the stable manifold is given by the surface y0+(1/3)(x0-y0)^2=0

% constants for the method
T = 79.1;
x0 = -1.11;
tol = 10^-7;
maxiter = 10000;

% array of num points
N = 100:150;

% is stable?
stable_jb = zeros(1,length(N));
stable_gs = zeros(1,length(N));
for j = 1:length(N);
    [~,~,~,m,ind] = fbeuler_jb( dx, dy, x0, T, N(j), tol, maxiter);
    p = polyfit(1:m,ind,1);
    if p(1) < 0;
        stable_jb(j) = 1;
    end
    [~,~,~,m,ind] = fbeuler_gs( dx, dy, x0, T, N(j), tol, maxiter);
    p = polyfit(1:m,ind,1);
    if p(1) < 0;
        stable_gs(j) = 1;
    end
    fprintf([num2str(j),'\n'])
end

plot(N,stable_jb,N,stable_gs)
axis([N(1),N(end),0,1.1])