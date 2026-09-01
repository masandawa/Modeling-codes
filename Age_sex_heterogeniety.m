%% ========================================================================
% STANDALONE AGE-SEX MALARIA HETEROGENEITY ANALYSIS
%
% Human groups:
%   FY = females aged <5 years
%   FA = females aged >=5 years
%   MY = males aged <5 years
%   MA = males aged >=5 years
%
% Human states:
%   S = susceptible
%   I = untreated infectious
%   H = formal treatment
%   M = self-medication
%   R = recovered
%
% Mosquito states:
%   Sv = susceptible mosquitoes
%   Ev = exposed mosquitoes
%   Iv = infectious mosquitoes
%
% Outputs:
%   1. Basic reproduction number R0
%   2. Population-wide endemic prevalence
%   3. Group-specific prevalence
%   4. Share of total prevalent infection burden
%   5. Human-to-mosquito infectious-pressure share
%   6. Female-male differences among >=5-year groups
%   7. Grouped-bar figure
%
% No external files are required.
% No results files are created.
% ========================================================================

clear;
clc;
close all;

%% ========================================================================
% REFERENCE PARAMETERS
% ========================================================================

p = reference_parameters();

%% ========================================================================
% DISEASE-FREE EQUILIBRIUM
% ========================================================================

[S0,~,~] = dfe_values(p);

%% ========================================================================
% INITIAL CONDITIONS
% ========================================================================

y0 = zeros(23,1);

% Human groups initially susceptible
for k = 1:4

    idx = 5*(k-1)+(1:5);

    y0(idx(1)) = S0(k);

end

% Small infection seed in each human group
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
% NUMERICAL SIMULATION
% ========================================================================

[t,Y] = ode15s( ...
    @(t,y) malaria_rhs(t,y,p), ...
    [0 Tend], ...
    y0, ...
    opts);

yEnd = Y(end,:)';

%% ========================================================================
% BASIC REPRODUCTION NUMBER
% ========================================================================

R0 = calc_R0(p);

%% ========================================================================
% ENDEMIC METRICS
% ========================================================================

out = calc_metrics(yEnd,p);

P    = out.P;
B    = out.B;
Pi   = out.Pi;
Pall = out.Pall;

%% ========================================================================
% FEMALE-MALE DIFFERENCES AMONG >=5-YEAR GROUPS
% ========================================================================

DeltaP_a = ...
    100*(P(2)-P(4));

DeltaPi_a = ...
    Pi(2)-Pi(4);

%% ========================================================================
% DISPLAY MAIN RESULTS
% ========================================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf('AGE-SEX HETEROGENEITY IN MALARIA BURDEN AND TRANSMISSION\n');
fprintf('=============================================================\n\n');

fprintf('Basic reproduction number, R0       = %.6f\n',R0);
fprintf('Population-wide endemic prevalence = %.6f %%\n',100*Pall);

fprintf('\n');
fprintf('Female-male prevalence difference among >=5 years:\n');
fprintf('Delta P_a = P_fa - P_ma = %.6f percentage points\n', ...
    DeltaP_a);

fprintf('\n');
fprintf('Female-male infectious-pressure difference among >=5 years:\n');
fprintf('Delta Pi_a = Pi_fa - Pi_ma = %.6f percentage points\n', ...
    DeltaPi_a);

%% ========================================================================
% RESULTS TABLE IN COMMAND WINDOW
% ========================================================================

PopulationGroup = { ...
    'Females <5 years'; ...
    'Females >=5 years'; ...
    'Males <5 years'; ...
    'Males >=5 years'};

Prevalence_percent = ...
    100*P;

BurdenShare_percent = ...
    B;

InfectiousPressureShare_percent = ...
    Pi;

ResultsTable = table( ...
    PopulationGroup, ...
    Prevalence_percent, ...
    BurdenShare_percent, ...
    InfectiousPressureShare_percent, ...
    'VariableNames', ...
    {'PopulationGroup', ...
     'Prevalence_percent', ...
     'BurdenShare_percent', ...
     'InfectiousPressureShare_percent'});

fprintf('\n');
fprintf('-------------------------------------------------------------\n');
fprintf('AGE-SEX-SPECIFIC RESULTS\n');
fprintf('-------------------------------------------------------------\n\n');

disp(ResultsTable);

%% ========================================================================
% ADDITIONAL SUMMARY
% ========================================================================

[~,idxPrev] = max(P);
[~,idxBurden] = max(B);
[~,idxPressure] = max(Pi);

shortNames = { ...
    'Female <5', ...
    'Female >=5', ...
    'Male <5', ...
    'Male >=5'};

fprintf('\n');
fprintf('Highest within-group prevalence      : %s\n', ...
    shortNames{idxPrev});

fprintf('Largest infection-burden share       : %s\n', ...
    shortNames{idxBurden});

fprintf('Largest infectious-pressure share    : %s\n', ...
    shortNames{idxPressure});

fprintf('\n');
fprintf('Combined burden share, >=5 years     = %.3f %%\n', ...
    B(2)+B(4));

fprintf('Combined pressure share, >=5 years   = %.3f %%\n', ...
    Pi(2)+Pi(4));

fprintf('\n');
fprintf('=============================================================\n');

%% ========================================================================
% GROUPED BAR FIGURE
% ========================================================================

plotData = [ ...
    100*P, ...
    B, ...
    Pi];

figure( ...
    'Color','w', ...
    'Position',[100 100 1200 650]);

bar(plotData,'grouped');

%% X-axis labels

ax = gca;

ax.XTick = 1:4;

ax.XTickLabel = { ...
    'Female $<5$', ...
    'Female $\geq5$', ...
    'Male $<5$', ...
    'Male $\geq5$'};

ax.TickLabelInterpreter = 'latex';

ax.FontSize = 11;
ax.LineWidth = 1;

%% Y-axis

ylabel( ...
    'Percentage (\%)', ...
    'Interpreter','latex', ...
    'FontSize',12);

ylim([0 60]);

%% Legend

legend( ...
    {'Within-group prevalence', ...
     'Share of total infection burden', ...
     'Share of infectious pressure'}, ...
    'Location','northwest', ...
    'FontSize',10);

%% Grid

grid on;
box on;

ax.XGrid = 'off';
ax.YGrid = 'on';
ax.GridAlpha = 0.15;

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

    % Group order:
    % 1 = FY
    % 2 = FA
    % 3 = MY
    % 4 = MA

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
% MALARIA ODE SYSTEM
% ========================================================================

function dy = malaria_rhs(~,y,p)

    dy = zeros(23,1);

    %% --------------------------------------------------------------
    % Human state matrix
    %
    % Rows:
    % 1 FY
    % 2 FA
    % 3 MY
    % 4 MA
    %
    % Columns:
    % 1 S
    % 2 I
    % 3 H
    % 4 M
    % 5 R
    % --------------------------------------------------------------

    V = reshape(y(1:20),5,4)';

    N = sum(V,2);

    %% Effective mosquito-contact factor

    q = ...
        p.chi(:) .* ...
        (1-p.epsN*p.c(:));

    Q = ...
        max(sum(q.*N),eps);

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

    %% Human infectious contribution

    infectious = ...
        sum( ...
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

    %% ==============================================================
    % HUMAN SYSTEM
    % ==============================================================

    for sex = 1:2

        if sex == 1

            ky = 1;
            ka = 2;

            Lam = p.Lambda_f;

        else

            ky = 3;
            ka = 4;

            Lam = p.Lambda_m;

        end

        %% ----------------------------------------------------------
        % <5-year states
        % ----------------------------------------------------------

        Sy = V(ky,1);
        Iy = V(ky,2);
        Hy = V(ky,3);
        My = V(ky,4);
        Ry = V(ky,5);

        %% ----------------------------------------------------------
        % >=5-year states
        % ----------------------------------------------------------

        Sa = V(ka,1);
        Ia = V(ka,2);
        Ha = V(ka,3);
        Ma = V(ka,4);
        Ra = V(ka,5);

        %% Treatment allocation

        tauHy = ...
            (1-p.pM(ky))*p.tauT;

        tauMy = ...
            p.pM(ky)*p.tauT;

        tauHa = ...
            (1-p.pM(ka))*p.tauT;

        tauMa = ...
            p.pM(ka)*p.tauT;

        %% ----------------------------------------------------------
        % <5-year equations
        % ----------------------------------------------------------

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

        %% ----------------------------------------------------------
        % >=5-year equations
        % ----------------------------------------------------------

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

        iy = ...
            5*(ky-1)+(1:5);

        ia = ...
            5*(ka-1)+(1:5);

        dy(iy) = dY;
        dy(ia) = dA;

    end

    %% ==============================================================
    % MOSQUITO SYSTEM
    % ==============================================================

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
% ENDEMIC METRICS
% ========================================================================

function out = calc_metrics(y,p)

    V = reshape(y(1:20),5,4)';

    N = sum(V,2);

    %% Effective mosquito-contact factor

    q = ...
        p.chi(:) .* ...
        (1-p.epsN*p.c(:));

    %% --------------------------------------------------------------
    % Prevalent infections
    %
    % D_g = I_g + H_g + M_g
    % --------------------------------------------------------------

    D = ...
        V(:,2) + ...
        V(:,3) + ...
        V(:,4);

    %% Within-group prevalence

    P = ...
        D ./ N;

    %% Share of total infection burden

    B = ...
        100*D/sum(D);

    %% --------------------------------------------------------------
    % Human-to-mosquito infectious pressure
    %
    % C_g = q_g(I_g + eta_H H_g + eta_M M_g)
    % --------------------------------------------------------------

    C = ...
        q .* ...
        ( ...
        V(:,2) + ...
        p.etaH*V(:,3) + ...
        p.etaM*V(:,4) ...
        );

    %% Share of total infectious pressure

    Pi = ...
        100*C/sum(C);

    %% Overall prevalence

    Pall = ...
        sum(D)/sum(N);

    %% Store

    out.N = N;
    out.D = D;

    out.P = P;
    out.B = B;

    out.C = C;
    out.Pi = Pi;

    out.Pall = Pall;

end

%% ========================================================================
% DISEASE-FREE EQUILIBRIUM
% ========================================================================

function [S0,Sv0,Q0] = dfe_values(p)

    S0 = zeros(4,1);

    %% Female <5 years

    S0(1) = ...
        p.Lambda_f / ...
        (p.mu_h+p.rho);

    %% Female >=5 years

    S0(2) = ...
        p.rho*S0(1) / ...
        p.mu_h;

    %% Male <5 years

    S0(3) = ...
        p.Lambda_m / ...
        (p.mu_h+p.rho);

    %% Male >=5 years

    S0(4) = ...
        p.rho*S0(3) / ...
        p.mu_h;

    %% Mosquito susceptible population

    Sv0 = ...
        p.Lambda_v / ...
        p.mu_v;

    %% Effective host availability

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

    %% ==============================================================
    % >=5-year infectious histories
    % ==============================================================

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

    %% ==============================================================
    % <5-year infectious histories
    % ==============================================================

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

    %% ==============================================================
    % COMMON HUMAN-MOSQUITO TRANSMISSION FACTOR
    % ==============================================================

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

    %% Contributions to R0^2

    contribution = ...
        factor * ...
        (q .* S0 .* T);

    %% Basic reproduction number

    R0 = ...
        sqrt(sum(contribution));

end