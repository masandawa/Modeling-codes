function S6_global_sensitivity_4panel
%% ========================================================================
% STANDALONE GLOBAL SENSITIVITY ANALYSIS
%
% Age-sex structured malaria model
%
% Human groups:
%   FY = females aged <5 years
%   FA = females aged >=5 years
%   MY = males aged <5 years
%   MA = males aged >=5 years
%
% Human states per group:
%   S, I, H, M, R
%
% Mosquito states:
%   Sv, Ev, Iv
%
% Global sensitivity:
%   Latin-hypercube sampling
%   Partial rank correlation coefficients (PRCC)
%   95% bootstrap confidence intervals
%
% Main manuscript outputs:
%   1. R0
%   2. Population-wide endemic prevalence
%   3. Delta P_a  = P_fa - P_ma
%   4. Delta Pi_a = Pi_fa - Pi_ma
%
% Main figure:
%   One 2 x 2 figure
%   Six strongest absolute PRCCs shown for each outcome
%
% This file is fully standalone.
% ========================================================================

clc;
close all;

%% ========================================================================
% REFERENCE PARAMETERS
% ========================================================================

p = reference_parameters();

%% ========================================================================
% INITIAL CONDITIONS
% ========================================================================

[S0,~,~] = dfe_values(p);

y0 = zeros(23,1);

% Human populations
for k = 1:4

    idx = 5*(k-1)+(1:5);

    y0(idx(1)) = S0(k);

end

% Small initial infection in each human group
seedH = 1e-4;

for k = 1:4

    idx = 5*(k-1)+(1:5);

    y0(idx(1)) = y0(idx(1)) - seedH;
    y0(idx(2)) = seedH;

end

% Mosquito initial conditions
y0(21) = p.m_v - 0.02;
y0(22) = 0.01;
y0(23) = 0.01;

%% ========================================================================
% NUMERICAL SETTINGS
% ========================================================================

Tend = 3650;

opts = odeset( ...
    'RelTol',1e-8, ...
    'AbsTol',1e-10, ...
    'NonNegative',1:23);

%% ========================================================================
% GLOBAL SENSITIVITY SETTINGS
% ========================================================================

rng(20260827);

N         = 1000;
rangeFrac = 0.20;
nBoot     = 500;

%% ========================================================================
% PARAMETERS INCLUDED IN GLOBAL SENSITIVITY
% ========================================================================

sensNames = { ...
    'a', ...
    'beta_vh', ...
    'beta_hv', ...
    'mu_v', ...
    'sigma_v', ...
    'm_v', ...
    'chi_fa', ...
    'chi_ma', ...
    'c_fa', ...
    'c_ma', ...
    'tauT', ...
    'pM_fa', ...
    'pM_ma', ...
    'gammaI', ...
    'gammaH', ...
    'gammaM', ...
    'etaH', ...
    'etaM'};

%% Mathematical labels for plotting

sensLabels = { ...
    '$a$', ...
    '$\beta_{vh}$', ...
    '$\beta_{hv}$', ...
    '$\mu_v$', ...
    '$\sigma_v$', ...
    '$m_v$', ...
    '$\chi_{fa}$', ...
    '$\chi_{ma}$', ...
    '$c_{fa}$', ...
    '$c_{ma}$', ...
    '$\tau_T$', ...
    '$p_{M,fa}$', ...
    '$p_{M,ma}$', ...
    '$\gamma_I$', ...
    '$\gamma_H$', ...
    '$\gamma_M$', ...
    '$\eta_H$', ...
    '$\eta_M$'};

%% Reference values

ref = [ ...
    p.a, ...
    p.beta_vh, ...
    p.beta_hv, ...
    p.mu_v, ...
    p.sigma_v, ...
    p.m_v, ...
    p.chi(2), ...
    p.chi(4), ...
    p.c(2), ...
    p.c(4), ...
    p.tauT, ...
    p.pM(2), ...
    p.pM(4), ...
    p.gammaI, ...
    p.gammaH, ...
    p.gammaM, ...
    p.etaH, ...
    p.etaM];

d = numel(sensNames);

%% ========================================================================
% LATIN-HYPERCUBE SAMPLING
%
% A local implementation is used so lhsdesign is not required.
% ========================================================================

U = lhs_local(N,d);

X = zeros(N,d);

for j = 1:d

    lo = (1-rangeFrac)*ref(j);
    hi = (1+rangeFrac)*ref(j);

    % Parameters naturally bounded between 0 and 1
    if ismember( ...
            sensNames{j}, ...
            {'c_fa','c_ma','pM_fa','pM_ma','etaH','etaM'})

        lo = max(0,lo);
        hi = min(1,hi);

    end

    X(:,j) = lo + (hi-lo)*U(:,j);

end

%% ========================================================================
% MODEL OUTPUTS
%
% Column 1 = R0
% Column 2 = overall prevalence (%)
% Column 3 = Delta P_a = P_fa - P_ma (percentage points)
% Column 4 = Delta Pi_a = Pi_fa - Pi_ma (percentage points)
% ========================================================================

Y = nan(N,4);

fprintf('\n=====================================================\n');
fprintf('GLOBAL SENSITIVITY ANALYSIS\n');
fprintf('=====================================================\n');
fprintf('Running %d sensitivity simulations...\n\n',N);

%% ========================================================================
% RUN SENSITIVITY SIMULATIONS
% ========================================================================

for i = 1:N

    ps = apply_sensitivity_row( ...
        p, ...
        X(i,:), ...
        sensNames);

    %% R0

    R0tmp = calc_R0(ps);

    %% Endemic equilibrium approximation

    oi = run_one( ...
        ps, ...
        y0, ...
        Tend, ...
        opts);

    %% Store outputs

    Y(i,1) = R0tmp;

    Y(i,2) = ...
        100*oi.Pall;

    Y(i,3) = ...
        100*(oi.P(2)-oi.P(4));

    Y(i,4) = ...
        oi.Pi(2)-oi.Pi(4);

    if mod(i,100)==0

        fprintf( ...
            'Completed %d of %d sensitivity simulations\n', ...
            i,N);

    end

end

%% ========================================================================
% OUTPUT NAMES
% ========================================================================

outputNames = { ...
    'R0', ...
    'OverallPrev_percent', ...
    'DeltaP_a_pp', ...
    'DeltaPi_a_pp'};

%% ========================================================================
% RANK TRANSFORMATION
%
% A local tied-rank function is used.
% ========================================================================

XR = zeros(size(X));
YR = zeros(size(Y));

for j = 1:d

    XR(:,j) = tiedrank_local(X(:,j));

end

for k = 1:4

    YR(:,k) = tiedrank_local(Y(:,k));

end

%% ========================================================================
% PRCC
% ========================================================================

PRCC = zeros(d,4);
PVAL = zeros(d,4);

for k = 1:4

    [PRCC(:,k),PVAL(:,k)] = ...
        prcc_one_output( ...
        XR, ...
        YR(:,k));

end

%% ========================================================================
% BOOTSTRAP CONFIDENCE INTERVALS
% ========================================================================

bootPRCC = nan(nBoot,d,4);

fprintf('\nRunning %d bootstrap replicates...\n\n',nBoot);

for b = 1:nBoot

    ii = randi(N,N,1);

    xb = XR(ii,:);

    for k = 1:4

        yb = YR(ii,k);

        bootPRCC(b,:,k) = ...
            prcc_one_output(xb,yb)';

    end

    if mod(b,100)==0

        fprintf( ...
            'Completed %d of %d bootstrap replicates\n', ...
            b,nBoot);

    end

end

%% ========================================================================
% 95% BOOTSTRAP CONFIDENCE INTERVALS
%
% Local percentile routine used to avoid toolbox dependence.
% ========================================================================

CIlo = zeros(d,4);
CIhi = zeros(d,4);

for k = 1:4

    for j = 1:d

        vals = squeeze(bootPRCC(:,j,k));

        CIlo(j,k) = percentile_local(vals,2.5);
        CIhi(j,k) = percentile_local(vals,97.5);

    end

end

%% ========================================================================
% DISPLAY COMPLETE PRCC RESULTS
% ========================================================================

for k = 1:4

    Tout = table( ...
        sensNames', ...
        PRCC(:,k), ...
        CIlo(:,k), ...
        CIhi(:,k), ...
        PVAL(:,k), ...
        'VariableNames', ...
        {'Parameter','PRCC','CI_low','CI_high','p_value'});

    [~,ord] = sort( ...
        abs(Tout.PRCC), ...
        'descend');

    Tout = Tout(ord,:);

    fprintf('\n=====================================================\n');
    fprintf('PRCC FOR %s\n',outputNames{k});
    fprintf('=====================================================\n');

    disp(Tout);

    writetable( ...
        Tout, ...
        ['PRCC_' outputNames{k} '.csv']);

end

%% ========================================================================
% IDENTIFY SIX STRONGEST PARAMETERS FOR EACH OUTPUT
% ========================================================================

nShow = 6;

ordTop = zeros(nShow,4);

for k = 1:4

    [~,ord] = sort( ...
        abs(PRCC(:,k)), ...
        'descend');

    ordTop(:,k) = ord(1:nShow);

end

%% ========================================================================
% CREATE COMPACT TABLE MATCHING THE FIGURE
% ========================================================================

Outcome   = cell(4*nShow,1);
Parameter = cell(4*nShow,1);

PRCCvalue = zeros(4*nShow,1);
CIlow     = zeros(4*nShow,1);
CIhigh    = zeros(4*nShow,1);
pvalue    = zeros(4*nShow,1);

outcomeTableNames = { ...
    'R0', ...
    'Population-wide prevalence', ...
    'DeltaP_a', ...
    'DeltaPi_a'};

row = 0;

for k = 1:4

    for j = 1:nShow

        row = row + 1;

        idx = ordTop(j,k);

        Outcome{row} = outcomeTableNames{k};
        Parameter{row} = sensNames{idx};

        PRCCvalue(row) = PRCC(idx,k);
        CIlow(row)     = CIlo(idx,k);
        CIhigh(row)    = CIhi(idx,k);
        pvalue(row)    = PVAL(idx,k);

    end

end

Tmain = table( ...
    Outcome, ...
    Parameter, ...
    PRCCvalue, ...
    CIlow, ...
    CIhigh, ...
    pvalue, ...
    'VariableNames', ...
    {'Outcome','Parameter','PRCC','CI_low','CI_high','p_value'});

fprintf('\n=====================================================\n');
fprintf('COMPACT TABLE USED IN THE MANUSCRIPT\n');
fprintf('=====================================================\n');

disp(Tmain);

writetable( ...
    Tmain, ...
    'PRCC_Main_Manuscript_Table.csv');

%% ========================================================================
% PUBLICATION-QUALITY FOUR-PANEL FIGURE
% ========================================================================

fig = figure( ...
    'Color','w', ...
    'Position',[80 50 1250 850]);

%% ------------------------------------------------------------------------
% PANEL (a): R0
% -------------------------------------------------------------------------

ax1 = subplot(2,2,1);

idx = ordTop(:,1);

barh( ...
    ax1, ...
    1:nShow, ...
    PRCC(idx,1), ...
    0.68);

hold(ax1,'on');

plot( ...
    ax1, ...
    [0 0], ...
    [0.5 nShow+0.5], ...
    'k-', ...
    'LineWidth',0.9);

hold(ax1,'off');

set(ax1, ...
    'YTick',1:nShow, ...
    'YTickLabel',sensLabels(idx), ...
    'TickLabelInterpreter','latex', ...
    'YDir','reverse', ...
    'FontSize',11, ...
    'LineWidth',1);

xlim(ax1,[-1 1]);
ylim(ax1,[0.5 nShow+0.5]);

xlabel(ax1,'PRCC','FontSize',11);

ylabel( ...
    ax1, ...
    'Parameter', ...
    'Interpreter','latex', ...
    'FontSize',11);

title( ...
    ax1, ...
    '(a) $\mathcal{R}_0$', ...
    'Interpreter','latex', ...
    'FontSize',12, ...
    'FontWeight','normal');

grid(ax1,'on');
box(ax1,'on');

ax1.XGrid = 'on';
ax1.YGrid = 'off';
ax1.GridAlpha = 0.15;

%% ------------------------------------------------------------------------
% PANEL (b): Overall prevalence
% -------------------------------------------------------------------------

ax2 = subplot(2,2,2);

idx = ordTop(:,2);

barh( ...
    ax2, ...
    1:nShow, ...
    PRCC(idx,2), ...
    0.68);

hold(ax2,'on');

plot( ...
    ax2, ...
    [0 0], ...
    [0.5 nShow+0.5], ...
    'k-', ...
    'LineWidth',0.9);

hold(ax2,'off');

set(ax2, ...
    'YTick',1:nShow, ...
    'YTickLabel',sensLabels(idx), ...
    'TickLabelInterpreter','latex', ...
    'YDir','reverse', ...
    'FontSize',11, ...
    'LineWidth',1);

xlim(ax2,[-1 1]);
ylim(ax2,[0.5 nShow+0.5]);

xlabel(ax2,'PRCC','FontSize',11);

ylabel( ...
    ax2, ...
    'Parameter', ...
    'Interpreter','latex', ...
    'FontSize',11);

title( ...
    ax2, ...
    '(b) Population-wide prevalence', ...
    'Interpreter','latex', ...
    'FontSize',12, ...
    'FontWeight','normal');

grid(ax2,'on');
box(ax2,'on');

ax2.XGrid = 'on';
ax2.YGrid = 'off';
ax2.GridAlpha = 0.15;

%% ------------------------------------------------------------------------
% PANEL (c): Delta P_a
% -------------------------------------------------------------------------

ax3 = subplot(2,2,3);

idx = ordTop(:,3);

barh( ...
    ax3, ...
    1:nShow, ...
    PRCC(idx,3), ...
    0.68);

hold(ax3,'on');

plot( ...
    ax3, ...
    [0 0], ...
    [0.5 nShow+0.5], ...
    'k-', ...
    'LineWidth',0.9);

hold(ax3,'off');

set(ax3, ...
    'YTick',1:nShow, ...
    'YTickLabel',sensLabels(idx), ...
    'TickLabelInterpreter','latex', ...
    'YDir','reverse', ...
    'FontSize',11, ...
    'LineWidth',1);

xlim(ax3,[-1 1]);
ylim(ax3,[0.5 nShow+0.5]);

xlabel(ax3,'PRCC','FontSize',11);

ylabel( ...
    ax3, ...
    'Parameter', ...
    'Interpreter','latex', ...
    'FontSize',11);

title( ...
    ax3, ...
    '(c) $\Delta P_a=P_{fa}-P_{ma}$', ...
    'Interpreter','latex', ...
    'FontSize',12, ...
    'FontWeight','normal');

grid(ax3,'on');
box(ax3,'on');

ax3.XGrid = 'on';
ax3.YGrid = 'off';
ax3.GridAlpha = 0.15;

%% ------------------------------------------------------------------------
% PANEL (d): Delta Pi_a
% -------------------------------------------------------------------------

ax4 = subplot(2,2,4);

idx = ordTop(:,4);

barh( ...
    ax4, ...
    1:nShow, ...
    PRCC(idx,4), ...
    0.68);

hold(ax4,'on');

plot( ...
    ax4, ...
    [0 0], ...
    [0.5 nShow+0.5], ...
    'k-', ...
    'LineWidth',0.9);

hold(ax4,'off');

set(ax4, ...
    'YTick',1:nShow, ...
    'YTickLabel',sensLabels(idx), ...
    'TickLabelInterpreter','latex', ...
    'YDir','reverse', ...
    'FontSize',11, ...
    'LineWidth',1);

xlim(ax4,[-1 1]);
ylim(ax4,[0.5 nShow+0.5]);

xlabel(ax4,'PRCC','FontSize',11);

ylabel( ...
    ax4, ...
    'Parameter', ...
    'Interpreter','latex', ...
    'FontSize',11);

title( ...
    ax4, ...
    '(d) $\Delta\Pi_a=\Pi_{fa}-\Pi_{ma}$', ...
    'Interpreter','latex', ...
    'FontSize',12, ...
    'FontWeight','normal');

grid(ax4,'on');
box(ax4,'on');

ax4.XGrid = 'on';
ax4.YGrid = 'off';
ax4.GridAlpha = 0.15;

%% ========================================================================
% EXPORT FIGURE
% ========================================================================

try

    exportgraphics( ...
        fig, ...
        'Figure_Global_Sensitivity_4Panel.png', ...
        'Resolution',300);

catch

    % Compatibility fallback for older MATLAB releases
    print( ...
        fig, ...
        'Figure_Global_Sensitivity_4Panel.png', ...
        '-dpng', ...
        '-r300');

end

%% ========================================================================
% SAVE COMPLETE NUMERICAL OUTPUT
% ========================================================================

save( ...
    'Global_Sensitivity_Results.mat', ...
    'X', ...
    'Y', ...
    'PRCC', ...
    'PVAL', ...
    'CIlo', ...
    'CIhi', ...
    'sensNames', ...
    'sensLabels', ...
    'outputNames', ...
    'ordTop', ...
    'Tmain');

fprintf('\n=====================================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('=====================================================\n');

fprintf('\nGenerated files:\n');
fprintf('  Figure_Global_Sensitivity_4Panel.png\n');
fprintf('  PRCC_Main_Manuscript_Table.csv\n');
fprintf('  PRCC_R0.csv\n');
fprintf('  PRCC_OverallPrev_percent.csv\n');
fprintf('  PRCC_DeltaP_a_pp.csv\n');
fprintf('  PRCC_DeltaPi_a_pp.csv\n');
fprintf('  Global_Sensitivity_Results.mat\n\n');

end


%% ========================================================================
% REFERENCE PARAMETERS
% ========================================================================

function p = reference_parameters()

p.mu_h    = 4.0e-5;
p.mu_v    = 0.052;
p.rho     = 1/(5*365);

p.sigma_v = 0.10;

p.a       = 0.50;
p.beta_vh = 0.086;
p.beta_hv = 0.48;

p.m_v     = 10;

p.omega   = 0.0027;
p.epsN    = 0.50;

% Order:
% FY, FA, MY, MA

p.chi = [ ...
    0.90, ...
    1.15, ...
    0.90, ...
    1.00];

p.c = [ ...
    0.50, ...
    0.169, ...
    0.50, ...
    0.234];

p.pM = [ ...
    0.05, ...
    0.260, ...
    0.05, ...
    0.426];

p.tauT = 0.060;

p.gammaI = 0.020;
p.gammaH = 0.20;
p.gammaM = 0.10;

p.etaH = 0.10;
p.etaM = 0.50;

p.delta = [ ...
    1.23e-5, ...
    1.53e-5, ...
    3.27e-5, ...
    1.23e-5];

p.Lambda_f = p.mu_h/2;
p.Lambda_m = p.mu_h/2;

p.Lambda_v = ...
    p.mu_v*p.m_v;

end


%% ========================================================================
% SIMULATE ONE PARAMETER SET
% ========================================================================

function out = run_one(p,y0,Tend,opts)

[~,Y] = ode15s( ...
    @(t,y) malaria_rhs(t,y,p), ...
    [0 Tend], ...
    y0, ...
    opts);

y = Y(end,:)';

out = calc_metrics(y,p);

end


%% ========================================================================
% MALARIA MODEL
% ========================================================================

function dy = malaria_rhs(~,y,p)

dy = zeros(23,1);

%% Human states
% Rows:
% FY, FA, MY, MA
%
% Columns:
% S, I, H, M, R

V = reshape(y(1:20),5,4)';

N = sum(V,2);

%% Effective mosquito contact

q = ...
    p.chi(:) .* ...
    (1-p.epsN*p.c(:));

Q = max(sum(q.*N),eps);

%% Mosquito states

Sv = y(21);
Ev = y(22);
Iv = y(23);

%% Human force of infection

lambda_h = ...
    p.a * ...
    p.beta_vh * ...
    q * ...
    Iv / Q;

%% Human infectious pressure

infectious = sum( ...
    q .* ...
    ( ...
    V(:,2) + ...
    p.etaH*V(:,3) + ...
    p.etaM*V(:,4) ...
    ));

%% Mosquito force of infection

lambda_v = ...
    p.a * ...
    p.beta_hv * ...
    infectious / Q;

%% Sex-specific human systems

for sex = 1:2

    if sex == 1

        ky  = 1;
        ka  = 2;
        Lam = p.Lambda_f;

    else

        ky  = 3;
        ka  = 4;
        Lam = p.Lambda_m;

    end

    %% Young states

    Sy = V(ky,1);
    Iy = V(ky,2);
    Hy = V(ky,3);
    My = V(ky,4);
    Ry = V(ky,5);

    %% >=5-year states

    Sa = V(ka,1);
    Ia = V(ka,2);
    Ha = V(ka,3);
    Ma = V(ka,4);
    Ra = V(ka,5);

    %% Treatment rates

    tauHy = ...
        (1-p.pM(ky))*p.tauT;

    tauMy = ...
        p.pM(ky)*p.tauT;

    tauHa = ...
        (1-p.pM(ka))*p.tauT;

    tauMa = ...
        p.pM(ka)*p.tauT;

    %% Individuals aged <5 years

    dY = [ ...

        Lam + ...
        p.omega*Ry - ...
        lambda_h(ky)*Sy - ...
        (p.mu_h+p.rho)*Sy;

        lambda_h(ky)*Sy - ...
        ( ...
        p.mu_h + ...
        p.rho + ...
        p.delta(ky) + ...
        p.gammaI + ...
        tauHy + ...
        tauMy ...
        )*Iy;

        tauHy*Iy - ...
        (p.mu_h+p.rho+p.gammaH)*Hy;

        tauMy*Iy - ...
        (p.mu_h+p.rho+p.gammaM)*My;

        p.gammaI*Iy + ...
        p.gammaH*Hy + ...
        p.gammaM*My - ...
        (p.mu_h+p.rho+p.omega)*Ry ...
        ];

    %% Individuals aged >=5 years

    dA = [ ...

        p.rho*Sy + ...
        p.omega*Ra - ...
        lambda_h(ka)*Sa - ...
        p.mu_h*Sa;

        p.rho*Iy + ...
        lambda_h(ka)*Sa - ...
        ( ...
        p.mu_h + ...
        p.delta(ka) + ...
        p.gammaI + ...
        tauHa + ...
        tauMa ...
        )*Ia;

        p.rho*Hy + ...
        tauHa*Ia - ...
        (p.mu_h+p.gammaH)*Ha;

        p.rho*My + ...
        tauMa*Ia - ...
        (p.mu_h+p.gammaM)*Ma;

        p.rho*Ry + ...
        p.gammaI*Ia + ...
        p.gammaH*Ha + ...
        p.gammaM*Ma - ...
        (p.mu_h+p.omega)*Ra ...
        ];

    iy = 5*(ky-1)+(1:5);
    ia = 5*(ka-1)+(1:5);

    dy(iy) = dY;
    dy(ia) = dA;

end

%% Mosquito system

dy(21) = ...
    p.Lambda_v - ...
    lambda_v*Sv - ...
    p.mu_v*Sv;

dy(22) = ...
    lambda_v*Sv - ...
    (p.sigma_v+p.mu_v)*Ev;

dy(23) = ...
    p.sigma_v*Ev - ...
    p.mu_v*Iv;

end


%% ========================================================================
% MODEL METRICS
% ========================================================================

function out = calc_metrics(y,p)

V = reshape(y(1:20),5,4)';

N = sum(V,2);

q = ...
    p.chi(:) .* ...
    (1-p.epsN*p.c(:));

%% Prevalent infection

D = ...
    V(:,2) + ...
    V(:,3) + ...
    V(:,4);

P = D./N;

%% Infectious pressure

C = ...
    q .* ...
    ( ...
    V(:,2) + ...
    p.etaH*V(:,3) + ...
    p.etaM*V(:,4) ...
    );

Pi = ...
    100*C/sum(C);

%% Outputs

out.P = P;
out.Pi = Pi;

out.Pall = ...
    sum(D)/sum(N);

end


%% ========================================================================
% DISEASE-FREE EQUILIBRIUM
% ========================================================================

function [S0,Sv0,Q0] = dfe_values(p)

S0 = zeros(4,1);

%% Female

S0(1) = ...
    p.Lambda_f / ...
    (p.mu_h+p.rho);

S0(2) = ...
    p.rho*S0(1) / ...
    p.mu_h;

%% Male

S0(3) = ...
    p.Lambda_m / ...
    (p.mu_h+p.rho);

S0(4) = ...
    p.rho*S0(3) / ...
    p.mu_h;

%% Mosquito

Sv0 = ...
    p.Lambda_v / ...
    p.mu_v;

%% Contact weighting

q = ...
    p.chi(:) .* ...
    (1-p.epsN*p.c(:));

Q0 = ...
    sum(q.*S0);

end


%% ========================================================================
% BASIC REPRODUCTION NUMBER
% ========================================================================

function R0 = calc_R0(p)

[S0,Sv0,Q0] = ...
    dfe_values(p);

q = ...
    p.chi(:) .* ...
    (1-p.epsN*p.c(:));

T  = zeros(4,1);
TH = zeros(4,1);
TM = zeros(4,1);

%% ------------------------------------------------------------------------
% Individuals aged >=5 years
% -------------------------------------------------------------------------

for ka = [2 4]

    tauH = ...
        (1-p.pM(ka))*p.tauT;

    tauM = ...
        p.pM(ka)*p.tauT;

    A = ...
        p.mu_h + ...
        p.delta(ka) + ...
        p.gammaI + ...
        tauH + ...
        tauM;

    BH = ...
        p.mu_h + ...
        p.gammaH;

    BM = ...
        p.mu_h + ...
        p.gammaM;

    TH(ka) = ...
        q(ka)*p.etaH/BH;

    TM(ka) = ...
        q(ka)*p.etaM/BM;

    T(ka) = ...
        q(ka)/A + ...
        (tauH/A)*TH(ka) + ...
        (tauM/A)*TM(ka);

end

%% ------------------------------------------------------------------------
% Individuals aged <5 years
% -------------------------------------------------------------------------

pair = [ ...
    1 2; ...
    3 4];

for z = 1:2

    ky = pair(z,1);
    ka = pair(z,2);

    tauH = ...
        (1-p.pM(ky))*p.tauT;

    tauM = ...
        p.pM(ky)*p.tauT;

    A = ...
        p.mu_h + ...
        p.rho + ...
        p.delta(ky) + ...
        p.gammaI + ...
        tauH + ...
        tauM;

    BH = ...
        p.mu_h + ...
        p.rho + ...
        p.gammaH;

    BM = ...
        p.mu_h + ...
        p.rho + ...
        p.gammaM;

    TH(ky) = ...
        q(ky)*p.etaH/BH + ...
        (p.rho/BH)*TH(ka);

    TM(ky) = ...
        q(ky)*p.etaM/BM + ...
        (p.rho/BM)*TM(ka);

    T(ky) = ...
        q(ky)/A + ...
        (tauH/A)*TH(ky) + ...
        (tauM/A)*TM(ky) + ...
        (p.rho/A)*T(ka);

end

%% Common factor

factor = ...
    p.a^2 * ...
    p.beta_vh * ...
    p.beta_hv * ...
    Sv0 * ...
    p.sigma_v / ...
    ( ...
    p.mu_v * ...
    (p.sigma_v+p.mu_v) * ...
    Q0^2 ...
    );

contribution = ...
    factor * ...
    (q.*S0.*T);

R0 = ...
    sqrt(sum(contribution));

end


%% ========================================================================
% APPLY SAMPLED PARAMETER VALUES
% ========================================================================

function ps = apply_sensitivity_row(ps,x,names)

for j = 1:numel(names)

    switch names{j}

        case 'a'
            ps.a = x(j);

        case 'beta_vh'
            ps.beta_vh = x(j);

        case 'beta_hv'
            ps.beta_hv = x(j);

        case 'mu_v'
            ps.mu_v = x(j);

        case 'sigma_v'
            ps.sigma_v = x(j);

        case 'm_v'
            ps.m_v = x(j);

        case 'chi_fa'
            ps.chi(2) = x(j);

        case 'chi_ma'
            ps.chi(4) = x(j);

        case 'c_fa'
            ps.c(2) = x(j);

        case 'c_ma'
            ps.c(4) = x(j);

        case 'tauT'
            ps.tauT = x(j);

        case 'pM_fa'
            ps.pM(2) = x(j);

        case 'pM_ma'
            ps.pM(4) = x(j);

        case 'gammaI'
            ps.gammaI = x(j);

        case 'gammaH'
            ps.gammaH = x(j);

        case 'gammaM'
            ps.gammaM = x(j);

        case 'etaH'
            ps.etaH = x(j);

        case 'etaM'
            ps.etaM = x(j);

    end

end

%% Maintain mosquito abundance parameterisation

ps.Lambda_v = ...
    ps.mu_v*ps.m_v;

end


%% ========================================================================
% PRCC FOR ONE OUTPUT
% ========================================================================

function [r,pv] = prcc_one_output(XR,y)

n = size(XR,1);
d = size(XR,2);

r  = zeros(d,1);
pv = zeros(d,1);

for j = 1:d

    others = setdiff(1:d,j);

    Z = [ ...
        ones(n,1), ...
        XR(:,others)];

    %% Residual parameter rank

    bx = Z\XR(:,j);

    rx = ...
        XR(:,j) - ...
        Z*bx;

    %% Residual output rank

    by = Z\y;

    ry = ...
        y - ...
        Z*by;

    %% Correlation

    denominator = ...
        sqrt(sum(rx.^2)*sum(ry.^2));

    if denominator <= eps

        r(j) = 0;

    else

        r(j) = ...
            sum(rx.*ry)/denominator;

    end

    %% p-value from Student t distribution

    df = ...
        n-d-1;

    tstat = ...
        r(j)*sqrt( ...
        df/max(1-r(j)^2,eps));

    pv(j) = ...
        student_t_two_sided_p(tstat,df);

end

end


%% ========================================================================
% TWO-SIDED STUDENT-t P-VALUE
%
% Uses betainc; Statistics Toolbox is not required.
% ========================================================================

function p = student_t_two_sided_p(t,nu)

x = ...
    nu/(nu+t.^2);

p = ...
    betainc(x,nu/2,0.5);

p = ...
    min(max(p,0),1);

end


%% ========================================================================
% LOCAL LATIN-HYPERCUBE SAMPLING
% ========================================================================

function U = lhs_local(N,d)

U = zeros(N,d);

for j = 1:d

    edges = ...
        ((0:N-1)' + rand(N,1))/N;

    order = randperm(N);

    U(:,j) = ...
        edges(order);

end

end


%% ========================================================================
% LOCAL TIED-RANK FUNCTION
% ========================================================================

function r = tiedrank_local(x)

x = x(:);

n = numel(x);

[xs,ord] = sort(x);

rs = zeros(n,1);

i = 1;

while i <= n

    j = i;

    while j < n && xs(j+1)==xs(i)

        j = j+1;

    end

    rankMean = ...
        (i+j)/2;

    rs(i:j) = ...
        rankMean;

    i = j+1;

end

r = zeros(n,1);

r(ord) = rs;

end


%% ========================================================================
% LOCAL PERCENTILE FUNCTION
% ========================================================================

function q = percentile_local(x,p)

x = sort(x(:));

x = x(~isnan(x));

n = numel(x);

if n == 0

    q = NaN;
    return;

end

if n == 1

    q = x;
    return;

end

pos = ...
    1 + (n-1)*(p/100);

lo = floor(pos);
hi = ceil(pos);

if lo == hi

    q = x(lo);

else

    q = ...
        x(lo) + ...
        (pos-lo)*(x(hi)-x(lo));

end

end