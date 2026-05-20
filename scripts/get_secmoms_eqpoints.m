% =========================================================================
% SECOND MOMENTS and EQUILIBRIUM POINTS 
% Computing the info for the histograms
% The equilibrium points were calculated with V=1000
% =========================================================================

% Results will be saved in the following folder:
% NOTE: This folder must be created manually before running this script.
output_folder = 'SecMoms_EqPoints/';


% Parameter scan: alpha3 values and parameter types
alpha3Vet = [0.03,0.3];
na        = length(alpha3Vet);
typeVet   = {'FIT', 'MIN'};  % Loop over both parameter types
nTypes    = length(typeVet);


% Load fitted parameters
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
y0    = [yy data(1,4) data(1,3) data(1,2)];
% Reduced model initial conditions (9 species: R3 and R3* removed)
y0_12 = [y0(1:2) y0(4:7) y0(9:11)];


% System dimensions
n = 11;
nSecondMoments = n*(n+1)/2;
nSecondMoments_12 = (n-2)*(n-2+1)/2;

% Integration time span
tspan = [3 60000];
V     = 1000;



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
    for i=1:na
        
        % Set alpha3 and alpha2
        params(10) = alpha3Vet(i);
        params(9) = rap_alphas*alpha3Vet(i);
    
    
        %--- COMPLETE MODEL ---
    
        % Solve Deterministic system and get equilibrium points
        [tsol,ysol] = ode45(@(t,y) sistemODE(t,y,params), tspan, y0); 
        equilibriumPoints = ysol(end,:);
        eq_Points = equilibriumPoints * V; % sol det * V
    
        % Save equilibrium points to file
        fileEqPoint = strcat(output_folder,type,'_eq_Points_a3_',num2str(params(10)),'_CM.dat');
        save(fileEqPoint,'eq_Points','-ascii');
    
    
        % Theoretical second moments
        [M,B] = calcMat_M_B(params,equilibriumPoints);
        % Build the linear system for the steady-state second moments
        [bigM, bigB] = matSecMoms(M,B,n);
        % Solve for the steady-state second moments
        secMom = -inv(bigM)*bigB;
    
        % Save second moments
        fileSecMom = strcat(output_folder,type,'_secMom_a3_',num2str(params(10)),'_CM.dat');
        save(fileSecMom,'secMom','-ascii')
    
    
    
        %--- REDUCED MODEL ---
    
        % Solve Deterministic system to get equilibrium points
        [tsol_12,ysol_12] = ode45(@(t,y) sistemODE_12(t,y,params), tspan, y0_12);
        equilibriumPoints_12 = ysol_12(end,:);
        eq_Points_12 = equilibriumPoints_12 * V; % sol det * V
    
        % Save equilibrium points to file
        fileEqPoint_12 = strcat(output_folder,type,'_eq_Points_a3_',num2str(params(10)),'_RM.dat');
        save(fileEqPoint_12,'eq_Points_12','-ascii');
     
        
        % Theoretical second moments
        [M_12,B_12] = calcMat_M_B_12(params,equilibriumPoints_12);
        % Build the linear system for the steady-state second moments
        [bigM_12, bigB_12] = matSecMoms(M_12,B_12,n-2);
        % Solve for the steady-state second moments
        secMom_12 = -inv(bigM_12)*bigB_12;
    
        % Save second moments
        fileSecMom_12 = strcat(output_folder,type,'_secMom_a3_',num2str(params(10)),'_RM.dat');
        save(fileSecMom_12,'secMom_12','-ascii')
    
    end

end