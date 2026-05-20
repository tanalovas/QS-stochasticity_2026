
function [t,sol]=solveLangevinLinear(M,G,nSteps,dt,startPoint,n)

% Solves a linear Langevin equation using the Euler-Maruyama scheme:
%
%   d(csi) = M*csi*dt + G*dW
%
% where M is the drift matrix, G is the noise matrix, and dW is a
% vector of independent Wiener increments (Gaussian white noise).
%
% Inputs:
%   M          - (n x n) drift matrix
%   G          - (n x n) noise matrix (square root of the diffusion matrix B)
%   nSteps     - number of time steps
%   dt         - time step size
%   startPoint - (n x 1) initial condition vector
%   n          - number of state variables (system dimension)
%
% Outputs:
%   t   - (1 x nSteps) time vector
%   sol - (nSteps x n) matrix of the solution trajectory (one row per time step)

% Build the time vector
t = 0:dt:dt*(nSteps-1);

% Pre-allocate the solution matrix
sol=zeros(nSteps,n);

% Set the initial condition
sol(1,:) = startPoint';
csiOld   = startPoint;

% Euler-Maruyama time integration
for i=2:nSteps

    % Normal random increments and apply noise matrix
    noise = G*randn(n, 1);

    % Update state: deterministic drift term + stochastic diffusion term
    csiNew = csiOld + dt*M*csiOld + sqrt(dt)*noise;
    
    % Next
    csiOld = csiNew;
    sol(i,:) = csiNew';

end

end
