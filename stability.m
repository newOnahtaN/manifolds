clear; clc; close all
% test stability on a linear test system

% evolution rules
A = 3;
d = 1;
dx = @(t,x,y) -A*x+d*y;
dy = @(t,x,y) A*y+d*x;

% fix initial conditions and so on
x0 = 10;
T = 50;
tol = 10^-10;
maxiter = 100;

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

plot(N,stable_jb,N,stable_gs); hold on;
gstability = (T+1)*(A+d)/2; % gershgorin stablility
plot(ones(1,length(N))*gstability, linspace(0,1.1,length(N)), '--k'); hold off
axis([N(1),N(end),0,1.1])

