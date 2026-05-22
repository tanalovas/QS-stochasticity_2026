% =========================================================================
% FITTING OF MODEL PARAMETERS TO EXPERIMENTAL DATA
% Fits the model parameters by minimising the error between the
% deterministic solution and the experimental time-series data for
% C8, HC10, and HC8. Runs nAverages independent optimisations and
% saves all results to a .dat file.
% =========================================================================


% Load experimental data
% Columns: time, C8_mean, HC10_mean, HC8_mean, C8_min, C8_max,
%          HC10_min, HC10_max, HC8_min, HC8_max
data = load('experimentalData.dat');


% Plot of raw experimental data
figure(1),
plot(data(:,1),data(:,2),'*-b',... % C8
     data(:,1),data(:,3),'*-r',... % HC10
     data(:,1),data(:,4),'*-k')    % HC8
legend('C8','HC10','HC8')



% --- Interpolation and time extension ---

%Interpolate first points
startT = data(1,1);
endT   = data(end,1);

% Time grid
xq   = startT:0.5:endT;

% Interpolate all columns onto the new grid
vq2  = interp1(data(:,1), data(:,2), xq);
vq3  = interp1(data(:,1), data(:,3), xq);
vq4  = interp1(data(:,1), data(:,4), xq);
vq5  = interp1(data(:,1), data(:,5), xq);
vq6  = interp1(data(:,1), data(:,6), xq);
vq7  = interp1(data(:,1), data(:,7), xq);
vq8  = interp1(data(:,1), data(:,8), xq);
vq9  = interp1(data(:,1), data(:,9), xq);
vq10 = interp1(data(:,1), data(:,10), xq);

% Extend the time axis and all signals by nPoints steps at constant value
% (flat extrapolation beyond the last data point, for ODE to settle)
nPoints = 50;

xq   = [xq   xq(end)+0.5:0.5:xq(end)+nPoints*0.5];
vq2  = [vq2  ones(1,nPoints)*vq2(end)];
vq3  = [vq3  ones(1,nPoints)*vq3(end)];
vq4  = [vq4  ones(1,nPoints)*vq4(end)];
vq5  = [vq5  ones(1,nPoints)*vq5(end)];
vq6  = [vq6  ones(1,nPoints)*vq6(end)];
vq7  = [vq7  ones(1,nPoints)*vq7(end)];
vq8  = [vq8  ones(1,nPoints)*vq8(end)];
vq9  = [vq9  ones(1,nPoints)*vq9(end)];
vq10 = [vq10 ones(1,nPoints)*vq10(end)];

% New matrix
data=[xq' vq2' vq3' vq4',vq5',vq6',vq7',vq8',vq9',vq10'];



% --- Optimisation settings ---

% Number of independent optimisation runs (different random starting points)
nAverages = 50;
dataForTheFile = zeros(nAverages,15);

% flagError controls which error metric is used by the objective function:
%   0 = global (all three species)
%   1 = C8 only
%   2 = HC10 only
%   3 = HC8 only
flagError  = 0;

% flagMethod selects the optimisation algorithm:
%   0 = SIMPLEX
%   1 = SIMPSA
flagMethod = 1;

% flagPlot controls whether the objective function produces intermediate plots
flagPlot   = 0;



% --- Optimisation loop ---
tic

for a=1:nAverages

    % Starting parameter guess
    P_subset = load('parLong.dat');
    
    % Upper and lower bounds for the fitted parameters 
    lb = 0.001*P_subset;
    ub = 1000*P_subset;
    

    if(flagMethod)
        % SIMPSA
        [ParBest,FVAL,EXITFLAG,OUTPUT] = SIMPSA('objectiveFunction3',P_subset,lb,ub,[],data,flagError,flagPlot);
    else
        % SIMPLEX
        options     = foptions;
        options(1)  = 0;
        options(2)  = 1.5;
        options(3)  = 1.5;
        options(14) =10000000;
        [ParBest,options] = SIMPLEXL('objectiveFunction3',P_subset,options,[],data,flagError);
    end

    
    % Save result only if the optimiser converged (EXITFLAG ~= 0)
    if(EXITFLAG)
        dataToBeSaved       = [ParBest FVAL];
        dataForTheFile(a,:) = dataToBeSaved;
    end
    
end

toc


% --- Save results ---

% File name encodes the error metric, method, and number of runs
fileName = sprintf('parFit_error_%d_method_%d_n_%d.dat',flagError,flagMethod,a);
save(fileName,'dataForTheFile','-ascii');

