% =========================================================================
% HISTOGRAMS OF EACH CASE
% Computing the info for the histograms
% The equilibrium points were calculated with V=1000
% =========================================================================

% This script reads input data from the following folder:
input_folder = 'SecMoms_EqPoints/';
% Figures will be saved in the following folder:
figures_folder = 'figures/';


%% Input data: equilibrium points and second moments

% All cases: {type, alpha3} combinations
cases = {'FIT', 0.044719; 'FIT', 0.4; 'MIN', 0.044719; 'MIN', 0.4};
nCases = size(cases, 1);

% Load all equilibrium points and second moments into cell arrays
eqPoints_CM = cell(nCases, 1);
eqPoints_RM = cell(nCases, 1);
secMom_CM   = cell(nCases, 1);
secMom_RM   = cell(nCases, 1);

for c = 1:nCases

    type  = cases{c,1};
    alpha3 = cases{c,2};

    tagEq  = strcat(input_folder,type,'_eq_Points_a3_',num2str(alpha3));
    tagSe  = strcat(input_folder,type,'_secMom_a3_',num2str(alpha3));

    eqPoints_CM{c} = load(strcat(tagEq,'CM.dat'));
    eqPoints_RM{c} = load(strcat(tagEq,'RM.dat'));
    secMom_CM{c}   = load(strcat(tagSe,'CM.dat'));
    secMom_RM{c}   = load(strcat(tagSe,'RM.dat'));

end


%% Simulation settings
nAverages = 10;
dt     = 1.25e-4;
nSteps = 9.6e8;

% Folders where Langevin output files are stored, one per case
langFolders_CM = {'FIT_0.044719/nAv_CM';
                  'FIT_0.4/nAv_CM';
                  'MIN_0.044719/nAv_CM';
                  'MIN_0.4/nAv_CM'};

langFolders_RM = {'FIT_0.044719/nAv_RM';
                  'FIT_0.4/nAv_RM';
                  'MIN_0.044719/nAv_RM';
                  'MIN_0.4/nAv_RM'};

% Load Langevin simulation data for all cases
stData_CM_all = cell(nCases, nAverages);
stData_RM_all = cell(nCases, nAverages);

for c = 1:nCases
    for nAv = 1:nAverages
        stData_CM_all{c, nAv} = load(strcat(langFolders_CM{c},'/fileOut_nAV_',num2str(nAv), '_CM.dat'));
        stData_RM_all{c, nAv} = load(strcat(langFolders_RM{c},'/fileOut_nAV_',num2str(nAv), '_RM.dat'));
    end
end


%% Indices and labels
sq_index    = [1,12,22,31,39,46,52,57,61,64,66];
sq_index_12 = [1,10,18,25,31,36,40,43,45];
index_12    = [1,2,4,5,6,7,9,10,11];
labelSp     = {'R1','R2','R3','M','I','R1*','R2*','R3*','HC8','HC10','C8'};
labelSp_12  = {'R1','R2','M','I','R1*','R2*','HC8','HC10','C8'};

nBins_CM = 50;
nBins_RM = 90;
nSp_CM   = 11;
nSp_RM   = 9;


%% Main loop: compute mean histogram across averages | all cases

for c = 1:nCases

    type  = cases{c,1};
    alpha3 = cases{c,2};

    % --- COMPLETE MODEL ---
    mean_bin_values = zeros(nSp_CM, nBins_CM);
    all_edges       = zeros(nSp_CM, nBins_CM);

    for sp = 1:nSp_CM
        min_bin = zeros(1, nAverages);
        max_bin = zeros(1, nAverages);

        % Find global bin range across all averages
        figure('visible', 'off');
        for plt = 1:nAverages
            y = histogram(stData_CM_all{c, plt}(:, sp+1), nBins_CM, 'FaceAlpha', .5);
            min_bin(plt) = min(y.BinEdges);
            max_bin(plt) = max(y.BinEdges);
        end
        close

        % Rebin all averages using the same edges
        bins_vector = linspace(min(min_bin), max(max_bin), nBins_CM+1);
        values_all  = zeros(nAverages, nBins_CM);
        figure('visible', 'off');
        for plt = 1:nAverages
            y = histogram(stData_CM_all{c, plt}(:, sp+1), 'BinEdges', bins_vector, 'FaceAlpha', .5);
            values_all(plt,:) = y.Values;
        end
        close

        mean_bin_values(sp,:) = mean(values_all, 1);
        all_edges(sp,:)       = bins_vector(1:end-1);

    end

    % --- REDUCED MODEL ---
    mean_bin_values_12 = zeros(nSp_RM, nBins_RM);
    all_edges_12       = zeros(nSp_RM, nBins_RM);

    for sp = 1:nSp_RM
        min_bin = zeros(1, nAverages);
        max_bin = zeros(1, nAverages);

        % Find global bin range across all averages
        figure('visible', 'off');
        for plt = 1:nAverages
            y = histogram(stData_RM_all{c, plt}(:, sp+1), nBins_RM, 'FaceAlpha', .5);
            min_bin(plt) = min(y.BinEdges);
            max_bin(plt) = max(y.BinEdges);
        end
        close

        % Rebin all averages using the same edges
        bins_vector = linspace(min(min_bin), max(max_bin), nBins_RM+1);
        values_all  = zeros(nAverages, nBins_RM);
        figure('visible', 'off');
        for plt = 1:nAverages
            y = histogram(stData_RM_all{c, plt}(:, sp+1), 'BinEdges', bins_vector, 'FaceAlpha', .5);
            values_all(plt,:) = y.Values;
        end
        close

        mean_bin_values_12(sp,:) = mean(values_all, 1);
        all_edges_12(sp,:)       = bins_vector(1:end-1);

    end


    % --- Plot ---
    meanval_CM_tog = mean_bin_values(index_12,:);
    edges_CM_tog   = all_edges(index_12,:);
    secmoms_CM_tog = secMom_CM{c}(sq_index);
    secmoms_CM_tog = secmoms_CM_tog(index_12);

    meanval_RM = mean_bin_values_12;
    edges_RM   = all_edges_12;
    secmoms_RM = secMom_RM{c}(sq_index_12);

    binWidth_CM_tog = edges_CM_tog(:,2) - edges_CM_tog(:,1);
    binWidth_RM     = edges_RM(:,2) - edges_RM(:,1);

    figure('visible', 'on');
    for sp = 1:nSp_RM

        % Normalization (area = 1)
        norm_mean_CM_tog = meanval_CM_tog(sp,:) / (sum(meanval_CM_tog(sp,:)) * binWidth_CM_tog(sp));
        norm_mean_RM     = meanval_RM(sp,:)     / (sum(meanval_RM(sp,:))     * binWidth_RM(sp));

        % Theoretical Gaussian curves
        sigma2_CM = secmoms_CM_tog(sp);
        sigma2_RM = secmoms_RM(sp);
        y_teoricoCM = 1/sqrt(2*pi*sigma2_CM) * exp(-(edges_CM_tog(sp,:).^2 / sigma2_CM / 2));
        y_teoricoRM = 1/sqrt(2*pi*sigma2_RM) * exp(-(edges_RM(sp,:).^2 / sigma2_RM / 2));

        subplot(3, 3, sp)
        % Langevin
        bar(edges_CM_tog(sp,:), norm_mean_CM_tog, 'FaceAlpha', .5), hold on
        bar(edges_RM(sp,:), norm_mean_RM, 'FaceAlpha', .5), hold on
        % Theoretical
        plot(edges_CM_tog(sp,:), y_teoricoCM, '--b', 'LineWidth', 1), hold on
        plot(edges_RM(sp,:), y_teoricoRM, '--r', 'LineWidth', 1)
        title(labelSp_12(sp))
        legend('CM', 'RM', 'T_{CM}', 'T_{RM}')

    end

    sgtitle(strcat(type,', alpha=',num2str(alpha3),', dt=',num2str(dt),', nSteps=',num2str(nSteps),', nAv ',num2str(nAverages)))
    name = strcat(figures_folder,'histogram_',type,'_alpha3_',num2str(alpha3),'.pdf');
    print(name, '-dpdf','-bestfit')

end
