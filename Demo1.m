clear all; clc; close all;
% test on the problem that I made

% define equations
dx = @(t,x,y) -x+2*y+(x-y).^2;
dy = @(t,x,y) (x-y).^2+y;
% the stable manifold is given by the surface y0+(1/3)(x0-y0)^2=0

% constants for the method
T = 4.1;
N = 310000;
x0 = -0.01;
tol = 10^-7;
maxiter = 10000;

% call method
[ x,y,t,m ] = fbeuler_jb( dx, dy, x0, T, N, tol, maxiter);