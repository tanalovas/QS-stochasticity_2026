% =========================================================================
% PLOT PARAMETER SCAN 
% =========================================================================

% This script reads input data from the following folder:
input_folder = 'scan_results/';
% Figures will be saved in the following folder:
figures_folder = 'figures/';


% INFO:

% alpha3
alpha3Vet=[0.03,0.3];

% all Kin and Katt that were used:
k_ins = [0.0002,0.0006,0.001,0.003,0.008,0.03,0.06,0.1,0.3,0.5,0.9];
k_atts = [0.0004,0.0008,0.002,0.005,0.009,0.01,0.03,0.06,0.1,0.2,0.5,1,1.5,2.5,5,7,10,13];


index = [1,2,4,5,6,7,9,10,11];
labelSp = {'[R1]','[R2]','[M]','[I]','[R1*]','[R2*]','[HC8]','[HC10]','[C8]'};


%% Ratios Second Moments


% Pre-allocate cell arrays to hold results for each alpha3 value
secmoms_CM_all = cell(length(alpha3Vet), 1);
secmoms_RM_all = cell(length(alpha3Vet), 1);

n=11;
comb = length(k_ins)*length(k_atts);

for at = 1:length(alpha3Vet)

    % Initializing 
    secmoms_CM = zeros(comb,n);
    secmoms_RM = zeros(comb,n-2);
    kin = zeros(comb,1);
    katt = zeros(comb,1);
    
    for i = 1:length(k_ins)
    
        for k = 1:length(k_atts)
    
            file_to_load_CM = strcat(input_folder,"scan_secMom_a3_",num2str(alpha3Vet(at)),"_kin_",num2str(k_ins(i)),"_katt_",num2str(k_atts(k)),"_CM.dat");
            secMom_v_CM=load(file_to_load_CM);
    
            file_to_load_RM = strcat(input_folder,"scan_secMom_a3_",num2str(alpha3Vet(at)),"_kin_",num2str(k_ins(i)),"_katt_",num2str(k_atts(k)),"_RM.dat");
            secMom_v_RM=load(file_to_load_RM);
    
            kin((i-1)*length(k_atts) + k) = k_ins(i);
            katt((i-1)*length(k_atts) + k) = k_atts(k);
    
            secmoms_CM((i-1)*length(k_atts) + k, :) = secMom_v_CM;
            secmoms_RM((i-1)*length(k_atts) + k, :) = secMom_v_RM;
    
        end
    end

    % Store the completed matrices for this alpha3 value
    secmoms_CM_all{at} = secmoms_CM;
    secmoms_RM_all{at} = secmoms_RM;

end

% Ratio of second moments to check variability between CM and RM:
secMoms_ratios_a3_1 = secmoms_CM_all{1}(:,index)./secmoms_RM_all{1};
secMoms_ratios_a3_2 = secmoms_CM_all{2}(:,index)./secmoms_RM_all{2};



%% Plot: Figure for the two alpha3 in HC8 HC10 C8

n=9;
last3sp = (n-2):n;  % indices of the last 3: HC8 HC10 C8

figure(1);
fig1 = tiledlayout(2, 3);

% alpha3 = 0.03
for i = last3sp
    nexttile;
    scatter(kin, katt, 20, secMoms_ratios_a3_1(:,i), 'filled')
    set(gca, 'XScale', 'log', 'YScale', 'log');
    colorbar
    title(labelSp(i))
    xticks(unique(kin))
    yticks(unique(katt))
    % reduce tick label font size
    ax = gca;
    ax.XAxis.FontSize = 8;
    ax.YAxis.FontSize = 5;
    % optionally rotate ticks
    xtickangle(80)
    % Add row label to the first tile of the row only
    if i == last3sp(1)
        ylabel(strcat('\alpha_3 = ', num2str(alpha3Vet(1))), 'FontSize', 11,'Rotation',0)
    end
end

% alpha3 = 0.3
for i = last3sp
    nexttile;
    scatter(kin, katt, 20, secMoms_ratios_a3_2(:,i), 'filled')
    set(gca, 'XScale', 'log', 'YScale', 'log');
    colorbar
    title(labelSp(i))
    xticks(unique(kin))
    yticks(unique(katt))
    % reduce tick label font size
    ax = gca;
    ax.XAxis.FontSize = 8;
    ax.YAxis.FontSize = 5;
    % optionally rotate ticks
    xtickangle(80)
        % Add row label to the first tile of the row only
    if i == last3sp(1)
        ylabel(strcat('\alpha_3 = ', num2str(alpha3Vet(2))), 'FontSize', 11,'Rotation',0)
    end
end

sgtitle('SecMoms ratio intensity');
xlabel(fig1, 'Kin');
ylabel(fig1, 'Katt');

name = strcat(figures_folder,'second_moments_ratio_intensity.pdf');
print(name, '-dpdf','-bestfit')



%% If you want to see all species:
% Choose secMoms_ratios_a3_1 or secMoms_ratios_a3_2


n=9;
figure(2);
fig2 = tiledlayout(3,3);
for i = 1:n
    nexttile;
    scatter(kin,katt,20,secMoms_ratios_a3_1(:,i),'filled')
    set(gca, 'XScale', 'log', 'YScale', 'log');
    colorbar
    title(labelSp(i))
    xticks(unique(kin))
    yticks(unique(katt))
    % reduce tick label font size
    ax = gca;
    ax.XAxis.FontSize = 8;
    ax.YAxis.FontSize = 5;
    % optionally rotate ticks
    xtickangle(80)
end
sgtitle(strcat('SecMoms ratio intensity with \alpha_3 = ',num2str(alpha3Vet(1))));
xlabel(fig2,'Kin');
ylabel(fig2,'Katt');

name = strcat(figures_folder,'second_moments_ratio_intensity_alpha3_',num2str(alpha3Vet(1)),'.pdf');
print(name, '-dpdf','-bestfit')