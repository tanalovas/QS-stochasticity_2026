% =========================================================================
% PLOT DETERMINISTIC SYSTEM
% =========================================================================

% Figures will be saved in the following folder:
figures_folder = 'figures/';


% Initial info
n=11;

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


% Load experimental data and set initial conditions
data = load('experimentalData.dat');
yy = paramFit(7:14); 
% Full model initial conditions (11 species)
y0 = [yy data(1,4) data(1,3) data(1,2)]; % HC8,HC10,C8
tspan = [3 31000];


% Notes:
% Experimental data are ordered: 
%   1      2       3         4       5      6        7       8       9      10
% time C8_mean HC10_mean HC8_mean C8_min C8_max HC10_min HC10_max HC8_min HC8_max
C8_neg = data(:,2) - data(:,5);
C8_pos = data(:,6) - data(:,2);
HC10_neg = data(:,3) - data(:,7);
HC10_pos = data(:,8) - data(:,3);
HC8_neg = data(:,4) - data(:,9);
HC8_pos = data(:,10) - data(:,4);

%% Deterministic Solutions
% COMPLETE MODEL

% Solve Deterministic system
[tsol,ysol] = ode45(@(t,y) sistemODE(t,y,params), tspan, y0); 

% Notes:
% Deterministic results (ysol) are ordered:
% [R1] [R2] [R3] [M] [I] [R1*] [R2*] [R3*] [HC8] [HC10] [C8]

% Plot versus experimental data
figure(1)

subplot(3,1,1)
plot(tsol,ysol(:,9),'k-','LineWidth',1.2); hold on % HC8
errorbar(data(:,1),data(:,4),HC8_neg,HC8_pos,'o','LineWidth',0.6,'Color',[0 0.5 1]); hold on
xlim([data(1,1),data(end,1)]); hold on
ylabel('[HC8]');
legend('Deterministic solution', 'Experimental data','Location','southeast')

subplot(3,1,2)
plot(tsol,ysol(:,10),'k-','LineWidth',1.2); hold on % HC10
errorbar(data(:,1),data(:,3),HC10_neg,HC10_pos,'o','LineWidth',0.6,'Color',[0 0.5 1]); hold on
xlim([data(1,1),data(end,1)]); hold on
ylabel('[HC10]');

subplot(3,1,3)
plot(tsol,ysol(:,11),'k-','LineWidth',1.2); hold on % C8
errorbar(data(:,1),data(:,2),C8_neg,C8_pos,'o','LineWidth',0.6,'Color',[0 0.5 1]); hold on
xlim([data(1,1),data(end,1)]); hold on
ylabel('[C8]');
xlabel('Time')

name = strcat(figures_folder,'deterministic_solution_vs_experimental_data.pdf');
print(name, '-dpdf','-bestfit')