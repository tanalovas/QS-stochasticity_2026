% =========================================================================
% PARAMETER SCAN 
% Over inhibition and activation rates kin and katt
% Test for different alpha 3 (scaling alpha 2 proportionally)
% =========================================================================

% Results will be saved in the following folder:
% NOTE: These folders must be created manually before running this script.
output_folder = 'scan_results/';
% Figures will be saved in the following folder:
figures_folder = 'scan_results/figures/';
 

tic

% Control flags to select which sections to run
indexFullModel=1;
indexReducedModel=1;
indexConvergence=0;


% Initial info
paramFit  = load('parLong.dat');
params = zeros(10,1);

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


% Parameter scan: alpha3 values
alpha3Vet = [params(10),0.4];
na        = length(alpha3Vet);


% Parameter scan: kin and katt grids
k_ins  = [0.0002,0.0006,0.001,0.003,0.008,0.03,0.06,0.1,0.3,0.5,0.9];
k_atts = [0.0004,0.0008,0.002,0.005,0.009,0.01,0.03,0.06,0.1,0.2,0.5,1,1.5,2.5,5,7,10,13];

nkin  = length(k_ins);
nkatt = length(k_atts);


% Indices of the diagonal entries of the second-moment vector
% (i.e. variances of each species), used to extract standard deviations
indexSecMom    = [1,12,22,31,39,46,52,57,61,64,66];
indexSecMom_12 = [1,10,18,25,31,36,40,43,45];



%% COMPLETE MODEL

if(indexFullModel)
    sprintf('Complete Model starts')

    for i=1:na
        
        % alpha3
        params(10) = alpha3Vet(i);
        sprintf(strcat('alpha3=',num2str(params(10))))

        % alpha2
        params(9) = rap_alphas*alpha3Vet(i);

        % Scan over kin values
        for kin = 1:nkin

            params(7) = k_ins(kin);

            % Scan over katt values
            for katt = 1:nkatt

                params(8) = k_atts(katt);
                sprintf(strcat('loop ',num2str(kin*katt*i), ' out of ',num2str(nkin*nkatt*na)))


                % Solve Deterministic system and get equilibrium points
                [tsol,ysol] = ode45(@(t,y) sistemODE(t,y,params), tspan, y0); 
                equilibriumPoints = ysol(end,:);

                % Save equilibrium points to file
                fileEqPoint = strcat(output_folder,'scan_eqPoints_a3_',num2str(params(10)),'_kin_',num2str(params(7)),'_katt_',num2str(params(8)),'_CM.dat');
                save(fileEqPoint,'equilibriumPoints','-ascii');


                % Theoretical second moments
                [M,B] = calcMat_M_B(params,equilibriumPoints);
                % Build the linear system for the steady-state second moments
                [bigM, bigB] = matSecMoms(M,B,n);
                % Solve for the steady-state second moments
                secMom = -inv(bigM)*bigB;

                % Extract standard deviations (square root of diagonal second moments)
                save_secMoms = zeros(1,n);
                for specie = 1:n
                    save_secMoms(specie) = sqrt(secMom(indexSecMom(specie)));
                end
                fileSecMom = strcat(output_folder,'scan_secMom_a3_',num2str(params(10)),'_kin_',num2str(params(7)),'_katt_',num2str(params(8)),'_CM.dat');
                save(fileSecMom,'save_secMoms','-ascii');



                % Check convergence plot for the ODE solution
                if(indexConvergence)
                    fig1 = figure('visible', 'off');
                    for plt = 1:n
                        subplot(3,4,plt);  
                        plot(tsol,ysol(:,plt),'b','LineWidth',2)
                    end
                    set(fig1,'PaperPositionMode','auto','PaperSize',[20,16])
                    name = strcat(figures_folder,'scan_convergence_K_eqPoints_complete_tspan_',num2str(tspan(2)),'_a3_',num2str(params(10)),'_kin_',num2str(params(7)),'_katt_',num2str(params(8)),'.pdf');
                    sgtitle(strcat('Convergence eqPoints Complete model (tspan=',num2str(tspan(2)),', alpha3=',num2str(params(10)),', kin=',num2str(params(7)),', katt=',num2str(params(8)),')'),'FontSize',20);
                    print(name, '-dpdf','-bestfit')
                    close(fig1)
                end

            end

        end

    end

end

sprintf('End of Complete model test')



%% REDUCED MODEL

if(indexReducedModel)
    sprintf('Reduced Model starts')


    for i=1:na
        
        % alpha3
        params(10) = alpha3Vet(i);

        %a lpha2
        params(9) = rap_alphas*alpha3Vet(i);
        
    
        % Scan over kin values
        for kin = 1:nkin
    
            params(7) = k_ins(kin);    

            % Scan over katt values
            for katt = 1:nkatt
                
                params(8) = k_atts(katt);
                sprintf(strcat('loop ',num2str(kin*katt*i), ' out of ',num2str(nkin*nkatt*na)))
    
                % Solve Deterministic system and get equilibrium points
                [tsol_12,ysol_12] = ode45(@(t,y) sistemODE_12(t,y,params), tspan, y0_12); 
                equilibriumPoints_12 = ysol_12(end,:);
    
                % Save equilibrium points
                fileEqPoint_12 = strcat(output_folder,'scan_eqPoints_a3_',num2str(params(10)),'_kin_',num2str(params(7)),'_katt_',num2str(params(8)),'_RM.dat');
                save(fileEqPoint_12,'equilibriumPoints_12','-ascii');
    
    
                % Theoretical second moments
                [M_12,B_12] = calcMat_M_B_12(params,equilibriumPoints_12);
                % Build the linear system for the steady-state second moments
                [bigM_12, bigB_12] = matSecMoms(M_12,B_12,n-2);
                % Solve for the steady-state second moments
                secMom_12 = -inv(bigM_12)*bigB_12;

                % Extract standard deviations (square root of diagonal second moments)
                save_secMoms_12 = zeros(1,n-2);
                for specie=1:n-2
                    save_secMoms_12(specie) = sqrt(secMom_12(indexSecMom_12(specie)));
                end
                fileSecMom_12 = strcat(output_folder,'scan_secMom_a3_',num2str(params(10)),'_kin_',num2str(params(7)),'_katt_',num2str(params(8)),'_RM.dat');
                save(fileSecMom_12,'save_secMoms_12','-ascii');
    
        
                % Check convergence plot for the ODE solution
                if(indexConvergence)
                    fig1 = figure('visible', 'off');
                    for plt = 1:n-2
                        subplot(3,3,plt);  
                        plot(tsol_12,ysol_12(:,plt),'b','LineWidth',2)
                    end
                    set(fig1,'PaperPositionMode','auto','PaperSize',[20,16])
                    name = strcat(figures_folder,'scan_convergence_K_eqPoints_reduced_tspan_',num2str(tspan(2)),'_a3_',num2str(params(10)),'_kin_',num2str(params(7)),'_katt_',num2str(params(8)),'.pdf');
                    sgtitle(strcat('Convergence eqPoints Reduced model (tspan=',num2str(tspan(2)),', alpha3=',num2str(params(10)),', kin=',num2str(params(7)),', katt=',num2str(params(8)),')'),'FontSize',20);
                    print(name, '-dpdf','-bestfit')
                    close(fig1)
                end
                
            end
    
        end
    
    end

end

sprintf('End of Reduced model test')

toc

