function [tsol,ysol] = solveODEfit(params,data)

tspan = data(1,1):0.01:data(end,1);

% R10 = 0.10;
% R20 = 0.10;
% R30 = 0.10;
% M0  = 0.10;
% I0  = 0.10;
% R1star0 = 0.10;
% R2star0 = 0.10;
% R3star0 = 0.10;
HC80  = data(1,4);
HC100 = data(1,3);
C80   = data(1,2);
yy    = params(7:14);
y0    = [yy HC80 HC100 C80];


[tsol,ysol] = ode45(@(t,y) sistemFit(t,y,params(1:6)), tspan, y0);

