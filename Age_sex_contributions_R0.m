%% ========================================================================
%  STANDALONE FIGURE:
%  Age-sex contributions to malaria transmission potential
%
%  Uses the manuscript percentages directly so the output matches
%  Table~\ref{tab:r0-contribution}.
%  No numerical labels are printed above the bars.
% ========================================================================

clear;
clc;
close all;

%% ------------------------------------------------------------------------
% Reference values reported in the manuscript
% -------------------------------------------------------------------------

Gamma = [ ...
    0.295, ...     % Females <5 years
    10.313, ...    % Females >=5 years
    0.294, ...     % Males <5 years
    7.550 ...      % Males >=5 years
    ];

% Percentage contributions reported from the full-precision analysis
Omega = [ ...
    1.60, ...      % Females <5 years
    55.89, ...     % Females >=5 years
    1.60, ...      % Males <5 years
    40.92 ...      % Males >=5 years
    ];

R0 = 4.2956;

groups = { ...
    'Females <5 years', ...
    'Females >=5 years', ...
    'Males <5 years', ...
    'Males >=5 years'};

%% ------------------------------------------------------------------------
% Command-window output
% -------------------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('AGE-SEX CONTRIBUTIONS TO MALARIA TRANSMISSION POTENTIAL\n');
fprintf('============================================================\n\n');

fprintf('%-22s %-20s %-18s\n', ...
    'Population group', 'R0^2 contribution', 'Percentage (%)');

fprintf('------------------------------------------------------------\n');

for k = 1:4
    fprintf('%-22s %-20.3f %-18.2f\n', ...
        groups{k}, Gamma(k), Omega(k));
end

fprintf('------------------------------------------------------------\n');

fprintf('Overall R0                = %.4f\n', R0);
fprintf('Combined >=5 contribution = %.2f%%\n', Omega(2) + Omega(4));

fprintf('============================================================\n\n');

%% ------------------------------------------------------------------------
% Figure
% -------------------------------------------------------------------------

figure('Color','w', ...
       'Position',[100 100 1100 650]);

bar(Omega,0.65);

%% Y-axis label
ylabel('Contribution to $\mathcal{R}_0^2$ (\%)', ...
    'Interpreter','latex', ...
    'FontSize',13);

%% X-axis categories
xticks(1:4);

xticklabels({ ...
    'Female $<5$', ...
    'Female $\geq 5$', ...
    'Male $<5$', ...
    'Male $\geq 5$'});

set(gca, ...
    'TickLabelInterpreter','latex', ...
    'FontSize',12, ...
    'LineWidth',1.0, ...
    'Box','on');

%% Axis range
ylim([0 60]);
yticks(0:10:60);

%% Grid
grid on;
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

%% Title
title('Age--sex contributions to malaria transmission potential', ...
    'Interpreter','latex', ...
    'FontSize',14, ...
    'FontWeight','normal');

