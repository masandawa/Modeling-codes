%% ========================================================================
% STANDALONE BEHAVIOURAL PATHWAY DECOMPOSITION
%
% Age-sex malaria model:
% behavioural pathways underlying female-male heterogeneity
%
% Scenarios:
%
% S0   = homogeneous behaviour
% SE   = exposure heterogeneity only
% SP   = prevention heterogeneity only
% ST   = treatment heterogeneity only
% Sall = all behavioural pathways
%
% Outcomes:
%
% Delta P_a  = P_fa - P_ma
% Delta Pi_a = Pi_fa - Pi_ma
%
% Positive values:
%   larger outcome among females aged >=5 years
%
% Negative values:
%   larger outcome among males aged >=5 years
%
% No external files are required.
% No files are saved.
% ========================================================================

clear;
clc;
close all;

%% ========================================================================
% NUMERICAL RESULTS FROM THE BEHAVIOURAL PATHWAY ANALYSIS
% ========================================================================

% Scenario order:
%
% 1 = S0
% 2 = SE
% 3 = SP
% 4 = ST
% 5 = Sall

Scenario = { ...
    'S0'; ...
    'SE'; ...
    'SP'; ...
    'ST'; ...
    'Sall'};

%% ------------------------------------------------------------------------
% Basic reproduction number
% -------------------------------------------------------------------------

R0 = [ ...
    4.287; ...     % S0
    4.297; ...     % SE
    4.287; ...     % SP
    4.287; ...     % ST
    4.296];        % Sall

%% ------------------------------------------------------------------------
% Female-male prevalence difference among >=5-year groups
%
% Delta P_a = P_fa - P_ma
%
% Units: percentage points
% -------------------------------------------------------------------------

DeltaP_a = [ ...
   -0.0001; ...    % S0
    0.0369; ...    % SE
    0.0094; ...    % SP
   -0.1444; ...    % ST
   -0.0975];       % Sall

%% ------------------------------------------------------------------------
% Female-male infectious-pressure difference among >=5-year groups
%
% Delta Pi_a = Pi_fa - Pi_ma
%
% Units: percentage points
% -------------------------------------------------------------------------

DeltaPi_a = [ ...
   -0.012; ...     % S0
    7.013; ...     % SE
    1.807; ...     % SP
   -1.833; ...     % ST
    7.014];        % Sall

%% ========================================================================
% DISPLAY RESULTS IN COMMAND WINDOW
% ========================================================================

fprintf('\n');
fprintf('===============================================================\n');
fprintf('BEHAVIOURAL PATHWAY DECOMPOSITION\n');
fprintf('===============================================================\n\n');

fprintf('%-12s %-10s %-18s %-18s\n', ...
    'Scenario','R0','Delta P_a (pp)','Delta Pi_a (pp)');

fprintf('---------------------------------------------------------------\n');

for k = 1:5

    fprintf('%-12s %-10.3f %-18.4f %-18.3f\n', ...
        Scenario{k}, ...
        R0(k), ...
        DeltaP_a(k), ...
        DeltaPi_a(k));

end

fprintf('---------------------------------------------------------------\n');

%% ========================================================================
% ADDITIVITY CHECK
%
% Effects are measured relative to homogeneous behaviour S0.
% ========================================================================

effect_E_P = DeltaP_a(2) - DeltaP_a(1);
effect_P_P = DeltaP_a(3) - DeltaP_a(1);
effect_T_P = DeltaP_a(4) - DeltaP_a(1);
effect_all_P = DeltaP_a(5) - DeltaP_a(1);

nonadd_P = ...
    effect_all_P - ...
    (effect_E_P + effect_P_P + effect_T_P);

effect_E_Pi = DeltaPi_a(2) - DeltaPi_a(1);
effect_P_Pi = DeltaPi_a(3) - DeltaPi_a(1);
effect_T_Pi = DeltaPi_a(4) - DeltaPi_a(1);
effect_all_Pi = DeltaPi_a(5) - DeltaPi_a(1);

nonadd_Pi = ...
    effect_all_Pi - ...
    (effect_E_Pi + effect_P_Pi + effect_T_Pi);

fprintf('\n');
fprintf('ADDITIVITY CHECK\n');
fprintf('---------------------------------------------------------------\n');

fprintf('Prevalence non-additivity          = %.5f pp\n', ...
    nonadd_P);

fprintf('Infectious-pressure non-additivity = %.5f pp\n', ...
    nonadd_Pi);

%% ========================================================================
% INTERPRETIVE SUMMARY
% ========================================================================

fprintf('\n');
fprintf('INTERPRETIVE SUMMARY\n');
fprintf('---------------------------------------------------------------\n');

fprintf(['Homogeneous behaviour produces almost no female-male ', ...
         'difference.\n']);

fprintf(['Exposure heterogeneity produces the largest isolated ', ...
         'positive shift in infectious pressure toward females.\n']);

fprintf(['Prevention heterogeneity produces a smaller shift in ', ...
         'the same direction.\n']);

fprintf(['Treatment heterogeneity shifts both prevalence and ', ...
         'infectious pressure toward males.\n']);

fprintf(['Under all pathways, prevalence is slightly higher among ', ...
         'males, whereas infectious pressure is higher among females.\n']);

fprintf('\n');
fprintf('===============================================================\n');

%% ========================================================================
% FIGURE
% ========================================================================

figure( ...
    'Color','w', ...
    'Position',[100 100 1200 520]);

%% Mathematical scenario labels

scenarioLabels = { ...
    '$\mathcal{S}_0$', ...
    '$\mathcal{S}_E$', ...
    '$\mathcal{S}_P$', ...
    '$\mathcal{S}_T$', ...
    '$\mathcal{S}_{\rm all}$'};

%% ========================================================================
% PANEL (a)
% Female-male prevalence difference
% ========================================================================

subplot(1,2,1);

bar( ...
    1:5, ...
    DeltaP_a, ...
    0.65);

hold on;

yline( ...
    0, ...
    'k-', ...
    'LineWidth',0.8);

hold off;

ax1 = gca;

ax1.XTick = 1:5;
ax1.XTickLabel = scenarioLabels;
ax1.TickLabelInterpreter = 'latex';

ax1.FontSize = 11;
ax1.LineWidth = 1;

ylabel( ...
    '$P_{fa}-P_{ma}$ (percentage points)', ...
    'Interpreter','latex', ...
    'FontSize',12);

title( ...
    '(a) Prevalence difference', ...
    'FontSize',12, ...
    'FontWeight','bold');

% Suitable range for values -0.1444 to +0.0369
ylim([-0.16 0.06]);

grid on;
box on;

ax1.XGrid = 'off';
ax1.YGrid = 'on';
ax1.GridAlpha = 0.15;

%% ========================================================================
% PANEL (b)
% Female-male infectious-pressure difference
% ========================================================================

subplot(1,2,2);

bar( ...
    1:5, ...
    DeltaPi_a, ...
    0.65);

hold on;

yline( ...
    0, ...
    'k-', ...
    'LineWidth',0.8);

hold off;

ax2 = gca;

ax2.XTick = 1:5;
ax2.XTickLabel = scenarioLabels;
ax2.TickLabelInterpreter = 'latex';

ax2.FontSize = 11;
ax2.LineWidth = 1;

ylabel( ...
    '$\Pi_{fa}-\Pi_{ma}$ (percentage points)', ...
    'Interpreter','latex', ...
    'FontSize',12);

title( ...
    '(b) Infectious-pressure difference', ...
    'FontSize',12, ...
    'FontWeight','bold');

% Suitable range for values -1.833 to +7.014
ylim([-2.5 8]);

grid on;
box on;

ax2.XGrid = 'off';
ax2.YGrid = 'on';
ax2.GridAlpha = 0.15;