% =========================================================================
% LANGEVIN vs THEORY 
% Comparison between stochastic simulations and analytical predictions 
% for the fluctuations of the signalling molecules
% =========================================================================

% Results will be saved in the folder structure listed below.
% NOTE: These folders must be created manually before running this script.
%
% FIT_0.03/nAv_CM
% FIT_0.03/nAv_RM
% FIT_0.3/nAv_CM
% FIT_0.3/nAv_RM
% MIN_0.03/nAv_CM
% MIN_0.03/nAv_RM
% MIN_0.3/nAv_CM
% MIN_0.3/nAv_RM
%


tic


% Initial info
n = 11;
% Number of independent Langevin trajectories to average over
nAverages = 25; 

% Parameter scan: alpha3 values and parameter types
alpha3Vet = [0.03,0.3];
na        = length(alpha3Vet);
typeVet   = {'FIT', 'MIN'};  % Loop over both parameter types
nTypes    = length(typeVet);


% Loading parameters
paramFit = load('parLong.dat');
params   = zeros(10,1);

% Basic transcription rate (beta)
params(1) = 1/80*60;

% Molecule degradation rate (mu)
params(2) = 0.001;

% Primary association constant (kon)
params(3) = 0.1*60;

% Primary dissociation constant (koff)
params(4) = 10*60;

% Secondary association constant (kon)
params(5) =  paramFit(1);

% Secondary dissociation constant (koff)
params(6) =  paramFit(2);

% Inhibition rate (kin)
params(7) = paramFit(3);

% Activation rate (katt)
params(8) = paramFit(4);

% R2 production rate (alpha2)
params(9) = paramFit(5);

% R3 production rate (alpha3)
params(10) = paramFit(6);

% Ratio alpha2/alpha3
rap_alphas = paramFit(5)/paramFit(6);


% Load experimental data and set initial conditions
data  = load('experimentalData.dat');
yy    = paramFit(7:14);
% Full model initial conditions (11 species)
y0    = [yy data(1,4) data(1,3) data(1,2)]; % HC8,HC10,C8
% Reduced model initial conditions (9 species: R3 and R3* removed)
y0_12 = [y0(1:2) y0(4:7) y0(9:11)];


% Integration time span for the deterministic ODE solver
tspan = [3 60000];

% Params for Langevin
dt     = 1.25e-4;
nSteps = 9.6e8;
skip   = 1e5;

% Saving info for histograms
setupdata = table(dt,nSteps,skip);
save('setup_info','setupdata');


%% Loop over parameter types (FIT and MIN)
for t = 1:nTypes

    type = typeVet{t};
    % Flag: 1 for MIN, 0 for FIT
    MIN  = strcmp(type, 'MIN');

    % Original values for inhibition and activation rate:
    params(7) = paramFit(3); %kin
    params(8) = paramFit(4); %katt

    % Configuration only if type = MIN
    if (MIN)
        params(7) = 0.0002; %kin
        params(8) = 13;     %katt
    end


    
    % Loop of alpha3
    for i = 1:na
    
        % alpha 3
        params(10) = alpha3Vet(i);
        
        % alpha 2
        params(9) = rap_alphas*alpha3Vet(i);
    
    
        %--- COMPLETE MODEL ---
        
        % Solve Deterministic system to evaluate equilibrium points
        [tsol,ysol] = ode45(@(t,y) sistemODE(t,y,params), tspan, y0); 
        equilibriumPoints = ysol(end,:);
        
        % Theoretical second moments
        [M,B] = calcMat_M_B(params,equilibriumPoints);
        
        % Get Matrix G for the noise:
        G=chol(B)'; % B=G*G'
        
        
        % Solve Langevin Equation
        for nAv=1:nAverages
        
            [t_out,sol] = solveLangevinLinear(M,G,nSteps,dt,zeros(n,1),n);
            % Subsample trajectory before saving to reduce file size
            matrixSol = [t_out(1:skip:end)' sol(1:skip:end,:)];
            
            fileOut = strcat(type,'_',num2str(params(10)),'/nAv_CM/fileOut_nAV_',num2str(nAv),'_CM.dat');
            save(fileOut,'matrixSol','-ascii');
        
        end
        
        
        %--- REDUCED MODEL ---
        
        % Solve Deterministic system to get equilibrium points
        [tsol_12,ysol_12] = ode45(@(t,y) sistemODE_12(t,y,params), tspan, y0_12);
        equilibriumPoints_12 = ysol_12(end,:);
        
        % Theoretical second moments
        [M_12,B_12] = calcMat_M_B_12(params,equilibriumPoints_12);
        
        
        % Get Matrix G for the noise:
        G_12=chol(B_12)'; % B=G*G'
        
        
        % Solve Langevin Equation
        for nAv=1:nAverages
        
            [t_12,sol_12] = solveLangevinLineare(M_12,G_12,nSteps,dt,zeros(n-2,1),n-2);
            % Subsample trajectory before saving to reduce file size
            matrixSol_12 = [t_12(1:skip:end)' sol_12(1:skip:end,:)];
            
            fileOut_12 = strcat(type,'_',num2str(params(10)),'/nAv_RM/fileOut_nAV_',num2str(nAv),'_RM.dat');
            save(fileOut_12,'matrixSol_12','-ascii');
        
        end
        
        
        display(strcat('Complete analysis for alpha3 = ',num2str(alpha3Vet(i))))
    
    
    
    end

    fprintf('Type of parameters: %s \n', type)

end



toc

