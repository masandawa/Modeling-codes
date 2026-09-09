%% ========================================================================
% CODE 5 -- FINAL FAST ROBUSTNESS ANALYSIS
%
% Robustness of seasonal malaria-control and resistance timing optima
%
% Examines:
%   A. seasonal forcing amplitude xi_v
%   B. effective SMC exposure q
%   C. number of SMC rounds K
%   D. resistant-strain transmission fitness f_R
%
% FAST SEARCH:
%   Stage 1 : 28-day coarse phase grid
%   Stage 2 : 2-day refinement around both preliminary minima
%   Stage 3 : 1-day local verification
%
% IMPORTANT:
%   A resistance optimum is NOT reported when resistance is effectively
%   phase-insensitive. NaN is retained internally but is not shown in the
%   manuscript figure.
%
% FINAL REFERENCE:
%   gamma = 1/20 day^-1
%   xi_v  = 0.40
%   q     = 0.70
%   K     = 4
%   tau   = 30 days
%   f_R   = 0.95
%   eps_R = 0.50
%
% ========================================================================

clear;
clc;
close all;

tic;

%% ========================================================================
% 1. REFERENCE PARAMETERS
% ========================================================================

base.Nh = 1e5;

base.a  = 0.33;

base.bS = 0.20;
base.cS = 0.40;

base.m  = 2.0;

base.sigma_h = 1/12;
base.sigma_v = 1/10;

base.gammaS = 1/20;
base.gammaR = 1/20;

base.mu_h = 1/(65*365);
base.mu_v = 0.0714;

%% Seasonality

base.T = 365;

base.xi_v = 0.40;

base.phi_v = 0;

base.LambdaBar = ...
    base.m*base.Nh*base.mu_v;

%% Resistant fitness

base.fR = 0.95;

base.bR = ...
    base.fR*base.bS;

base.cR = ...
    base.fR*base.cS;

%% SMC programme

base.q = 0.70;

base.K = 4;

base.tau = 30;

%% Drug response

base.eS = 0.98;
base.eR = 0.90;

base.epsS = 0.90;
base.epsR = 0.50;

base.DP = 35;

base.omega = ...
    1/base.DP;

%% Resistance seeding

base.pR0 = 0.01;

%% Simulation horizon

base.years = 20;

base.evalYears = 5;

%% ========================================================================
% 2. FAST NUMERICAL SETTINGS
% ========================================================================

base.opts6 = odeset( ...
    'RelTol',1e-6, ...
    'AbsTol',1e-8, ...
    'MaxStep',5, ...
    'NonNegative',1:6);

base.opts11 = odeset( ...
    'RelTol',1e-6, ...
    'AbsTol',1e-8, ...
    'MaxStep',5, ...
    'NonNegative',1:11);

%% Periodic resident convergence

base.periodicTol = 1e-7;

base.maxBurnYears = 30;

%% Phase-identifiability thresholds

base.burdenFlatTolerancePP = 0.01;

base.resistanceFlatTolerancePP = 0.10;

%% ========================================================================
% 3. REPRODUCTION NUMBERS
% ========================================================================

R0S = reproduction_number( ...
    base,base.bS,base.cS,base.gammaS);

R0R = reproduction_number( ...
    base,base.bR,base.cR,base.gammaR);

fprintf('\n============================================================\n');
fprintf('CODE 5: FAST ROBUSTNESS OF SEASONAL TIMING OPTIMA\n');
fprintf('============================================================\n');

fprintf('Reference R0S = %.6f\n',R0S);
fprintf('Reference R0R = %.6f\n',R0R);

fprintf('Reference xi_v = %.2f\n',base.xi_v);
fprintf('Reference q    = %.2f\n',base.q);
fprintf('Reference K    = %d\n',base.K);
fprintf('Reference f_R  = %.2f\n',base.fR);

%% ========================================================================
% 4. ROBUSTNESS SCENARIOS
% ========================================================================

xiVals = ...
    [0.00 0.20 0.40 0.60];

qVals = ...
    [0.50 0.70 0.90];

KVals = ...
    [3 4 5];

fitVals = ...
    [0.90 0.95 1.00];

%% ========================================================================
% 5. STORAGE
% ========================================================================

A = initialise_results(length(xiVals));
B = initialise_results(length(qVals));
C = initialise_results(length(KVals));
D = initialise_results(length(fitVals));

%% ========================================================================
% 6. REFERENCE UNTREATED PERIODIC STATE
%
% Reused for q, K and f_R scenarios.
% ========================================================================

fprintf('\nComputing untreated reference periodic state...\n');

[zReference,nBurn,burnErr] = ...
    periodic_resident(base);

fprintf( ...
    'Periodic state converged after %d annual cycles; error = %.3e\n', ...
    nBurn,burnErr);

%% ========================================================================
% 7. A. TRANSMISSION SEASONALITY
% ========================================================================

fprintf('\n============================================================\n');
fprintf('A. TRANSMISSION SEASONALITY\n');
fprintf('============================================================\n');

for j = 1:length(xiVals)

    p = base;

    p.xi_v = ...
        xiVals(j);

    %% Non-seasonal case has no meaningful calendar-phase optimum

    if abs(p.xi_v) < 1e-12

        S.phiC = NaN;
        S.phiR = NaN;
        S.dCR  = NaN;

        S.burdenRangePP     = 0;
        S.resistanceRangePP = 0;

        S.controlPenalty       = NaN;
        S.resistancePenaltyPP  = NaN;

        S.controlFlat    = true;
        S.resistanceFlat = true;

    else

        if abs(p.xi_v-base.xi_v) < 1e-12

            z0 = ...
                zReference;

        else

            [z0,~,~] = ...
                periodic_resident(p);

        end

        S = ...
            fast_scenario_optima(p,z0);

    end

    A = ...
        store_result(A,j,S);

    print_result( ...
        sprintf('xi_v = %.2f',xiVals(j)), ...
        A,j);

end

%% ========================================================================
% 8. B. EFFECTIVE SMC EXPOSURE
% ========================================================================

fprintf('\n============================================================\n');
fprintf('B. EFFECTIVE SMC EXPOSURE\n');
fprintf('============================================================\n');

for j = 1:length(qVals)

    p = base;

    p.q = ...
        qVals(j);

    S = ...
        fast_scenario_optima( ...
        p,zReference);

    B = ...
        store_result(B,j,S);

    print_result( ...
        sprintf('q = %.2f',qVals(j)), ...
        B,j);

end

%% ========================================================================
% 9. C. NUMBER OF SMC ROUNDS
% ========================================================================

fprintf('\n============================================================\n');
fprintf('C. NUMBER OF SMC ROUNDS\n');
fprintf('============================================================\n');

for j = 1:length(KVals)

    p = base;

    p.K = ...
        KVals(j);

    S = ...
        fast_scenario_optima( ...
        p,zReference);

    C = ...
        store_result(C,j,S);

    print_result( ...
        sprintf('K = %d',KVals(j)), ...
        C,j);

end

%% ========================================================================
% 10. D. RESISTANT-STRAIN FITNESS
% ========================================================================

fprintf('\n============================================================\n');
fprintf('D. RESISTANT-STRAIN FITNESS\n');
fprintf('============================================================\n');

for j = 1:length(fitVals)

    p = base;

    p.fR = ...
        fitVals(j);

    p.bR = ...
        p.fR*p.bS;

    p.cR = ...
        p.fR*p.cS;

    S = ...
        fast_scenario_optima( ...
        p,zReference);

    D = ...
        store_result(D,j,S);

    print_result( ...
        sprintf('f_R = %.2f',fitVals(j)), ...
        D,j);

end

%% ========================================================================
% 11. RESULTS TABLES
% ========================================================================

SeasonalityResults = ...
    make_result_table( ...
    xiVals,A,'xi_v');

ExposureResults = ...
    make_result_table( ...
    qVals,B,'q');

RoundsResults = ...
    make_result_table( ...
    KVals,C,'K');

FitnessResults = ...
    make_result_table( ...
    fitVals,D,'f_R');

fprintf('\n============================================================\n');
fprintf('SEASONALITY RESULTS\n');
fprintf('============================================================\n');
disp(SeasonalityResults);

fprintf('\n============================================================\n');
fprintf('SMC EXPOSURE RESULTS\n');
fprintf('============================================================\n');
disp(ExposureResults);

fprintf('\n============================================================\n');
fprintf('NUMBER-OF-ROUNDS RESULTS\n');
fprintf('============================================================\n');
disp(RoundsResults);

fprintf('\n============================================================\n');
fprintf('RESISTANT-FITNESS RESULTS\n');
fprintf('============================================================\n');
disp(FitnessResults);

%% ========================================================================
% 12. CLEAN MANUSCRIPT FIGURE
%
% No NI text.
% No d_CR labels.
% No overlapping annotations.
%
% Missing phi_R* values simply appear as gaps.
% ========================================================================

figure( ...
    'Color','w', ...
    'Position',[70 60 1250 780]);

TL = tiledlayout( ...
    2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

%% ------------------------------------------------------------------------
% PANEL A
% -------------------------------------------------------------------------

ax1 = nexttile;

plot_optima( ...
    ax1, ...
    xiVals, ...
    A, ...
    'Seasonal forcing amplitude, \xi_v', ...
    '(a) Transmission seasonality');

xlim(ax1,[0.18 0.62]);

%% ------------------------------------------------------------------------
% PANEL B
% -------------------------------------------------------------------------

ax2 = nexttile;

plot_optima( ...
    ax2, ...
    qVals, ...
    B, ...
    'Effective SMC exposure, q', ...
    '(b) Intervention exposure');

xlim(ax2,[0.48 0.92]);

%% ------------------------------------------------------------------------
% PANEL C
% -------------------------------------------------------------------------

ax3 = nexttile;

plot_optima( ...
    ax3, ...
    KVals, ...
    C, ...
    'Number of SMC rounds, K', ...
    '(c) Treatment frequency');

xlim(ax3,[2.8 5.2]);

xticks(ax3,KVals);

%% ------------------------------------------------------------------------
% PANEL D
% -------------------------------------------------------------------------

ax4 = nexttile;

plot_optima( ...
    ax4, ...
    fitVals, ...
    D, ...
    'Relative resistant fitness, f_R', ...
    '(d) Resistant-strain fitness');

xlim(ax4,[0.89 1.01]);

%% ------------------------------------------------------------------------
% SAME AXIS SCALE IN ALL PANELS
% -------------------------------------------------------------------------

allAxes = ...
    [ax1 ax2 ax3 ax4];

for ax = allAxes

    ylim(ax,[20 150]);

    yticks(ax,20:20:140);

    grid(ax,'on');

    box(ax,'on');

    ax.FontSize = 10;

    ax.LineWidth = 0.8;

end

%% ------------------------------------------------------------------------
% ONE SHARED LEGEND
% -------------------------------------------------------------------------

lgd = legend( ...
    ax1, ...
    {'Malaria-control optimum, \phi_C^*', ...
     'Resistance optimum, \phi_R^*'}, ...
    'Orientation','horizontal', ...
    'FontSize',9);

lgd.Layout.Tile = 'south';

%% ------------------------------------------------------------------------
% OVERALL TITLE
% -------------------------------------------------------------------------



%% ========================================================================
% 13. RUNTIME
% ========================================================================

elapsedTime = toc;

fprintf('\n============================================================\n');
fprintf('TOTAL RUNTIME = %.2f seconds\n',elapsedTime);
fprintf('============================================================\n');

%% ========================================================================
% OPTIONAL EXPORT
% ========================================================================

% exportgraphics( ...
%     gcf, ...
%     'Timing_robustness_final.png', ...
%     'Resolution',600);

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================


%% ------------------------------------------------------------------------
% Reproduction number
% -------------------------------------------------------------------------

function R0 = reproduction_number(p,b,c,gamma)

R0 = sqrt( ...
    (p.m*p.a^2*b*c / ...
    ((gamma+p.mu_h)*p.mu_v)) ...
    *(p.sigma_h/(p.sigma_h+p.mu_h)) ...
    *(p.sigma_v/(p.sigma_v+p.mu_v)) );

end


%% ------------------------------------------------------------------------
% Initialize results
% -------------------------------------------------------------------------

function R = initialise_results(n)

R.phiC = nan(1,n);
R.phiR = nan(1,n);
R.dCR  = nan(1,n);

R.burdenRangePP     = nan(1,n);
R.resistanceRangePP = nan(1,n);

R.controlPenalty      = nan(1,n);
R.resistancePenaltyPP = nan(1,n);

R.controlFlat    = false(1,n);
R.resistanceFlat = false(1,n);

end


%% ------------------------------------------------------------------------
% Store results
% -------------------------------------------------------------------------

function R = store_result(R,j,S)

R.phiC(j) = S.phiC;
R.phiR(j) = S.phiR;
R.dCR(j)  = S.dCR;

R.burdenRangePP(j) = ...
    S.burdenRangePP;

R.resistanceRangePP(j) = ...
    S.resistanceRangePP;

R.controlPenalty(j) = ...
    S.controlPenalty;

R.resistancePenaltyPP(j) = ...
    S.resistancePenaltyPP;

R.controlFlat(j) = ...
    S.controlFlat;

R.resistanceFlat(j) = ...
    S.resistanceFlat;

end


%% ------------------------------------------------------------------------
% Console output
% -------------------------------------------------------------------------

function print_result(label,R,j)

fprintf('\n%s\n',label);

if R.controlFlat(j)

    fprintf('  malaria-control optimum: phase-insensitive\n');

else

    fprintf( ...
        '  phi_C* = %.0f d\n', ...
        R.phiC(j));

end

if R.resistanceFlat(j)

    fprintf('  resistance optimum: phase-insensitive\n');

else

    fprintf( ...
        '  phi_R* = %.0f d\n', ...
        R.phiR(j));

end

if ~isnan(R.dCR(j))

    fprintf( ...
        '  d_CR = %.1f d\n', ...
        R.dCR(j));

end

fprintf( ...
    '  malaria phase range    = %.3f pp\n', ...
    R.burdenRangePP(j));

fprintf( ...
    '  resistance phase range = %.3f pp\n', ...
    R.resistanceRangePP(j));

end


%% ------------------------------------------------------------------------
% Results table
% -------------------------------------------------------------------------

function T = make_result_table(values,R,varName)

T = table( ...
    values(:), ...
    R.phiC(:), ...
    R.phiR(:), ...
    R.dCR(:), ...
    R.burdenRangePP(:), ...
    R.resistanceRangePP(:), ...
    R.controlPenalty(:), ...
    R.resistancePenaltyPP(:), ...
    R.controlFlat(:), ...
    R.resistanceFlat(:), ...
    'VariableNames',{ ...
    varName, ...
    'phiC_day', ...
    'phiR_day', ...
    'dCR_day', ...
    'MalariaPhaseRange_pp', ...
    'ResistancePhaseRange_pp', ...
    'ControlPenalty_percent', ...
    'ResistancePenalty_pp', ...
    'ControlPhaseInsensitive', ...
    'ResistancePhaseInsensitive'});

end


%% ========================================================================
% FAST THREE-STAGE PHASE SEARCH
% ========================================================================

function S = fast_scenario_optima(p,z0)

%% ------------------------------------------------------------------------
% Stage 1: 28-day search
% -------------------------------------------------------------------------

phi1 = ...
    0:28:364;

P1 = ...
    nan(size(phi1));

R1 = ...
    nan(size(phi1));

for i = 1:length(phi1)

    [P1(i),R1(i)] = ...
        evaluate_phase( ...
        phi1(i),z0,p);

end

[~,iC1] = ...
    min(P1);

[~,iR1] = ...
    min(R1);

phiC1 = ...
    phi1(iC1);

phiR1 = ...
    phi1(iR1);

%% ------------------------------------------------------------------------
% Stage 2: 2-day refinement
% -------------------------------------------------------------------------

phi2 = ...
    unique( ...
    mod( ...
    [phiC1+(-16:2:16), ...
     phiR1+(-16:2:16)], ...
    p.T));

P2 = ...
    nan(size(phi2));

R2 = ...
    nan(size(phi2));

for i = 1:length(phi2)

    [P2(i),R2(i)] = ...
        evaluate_phase( ...
        phi2(i),z0,p);

end

[~,iC2] = ...
    min(P2);

[~,iR2] = ...
    min(R2);

phiC2 = ...
    phi2(iC2);

phiR2 = ...
    phi2(iR2);

%% ------------------------------------------------------------------------
% Stage 3: 1-day verification
%
% +/-3 days makes the final minimum less sensitive to the preceding
% 2-day grid while remaining computationally inexpensive.
% -------------------------------------------------------------------------

phi3 = ...
    unique( ...
    mod( ...
    [phiC2+(-3:3), ...
     phiR2+(-3:3)], ...
    p.T));

P3 = ...
    nan(size(phi3));

R3 = ...
    nan(size(phi3));

for i = 1:length(phi3)

    [P3(i),R3(i)] = ...
        evaluate_phase( ...
        phi3(i),z0,p);

end

%% ------------------------------------------------------------------------
% Overall phase ranges
% -------------------------------------------------------------------------

Pall = ...
    [P1 P2 P3];

Rall = ...
    [R1 R2 R3];

S.burdenRangePP = ...
    100*(max(Pall)-min(Pall));

S.resistanceRangePP = ...
    100*(max(Rall)-min(Rall));

%% ------------------------------------------------------------------------
% Phase sensitivity classification
% -------------------------------------------------------------------------

S.controlFlat = ...
    S.burdenRangePP < ...
    p.burdenFlatTolerancePP;

S.resistanceFlat = ...
    S.resistanceRangePP < ...
    p.resistanceFlatTolerancePP;

%% ------------------------------------------------------------------------
% Final control optimum
% -------------------------------------------------------------------------

if S.controlFlat

    S.phiC = ...
        NaN;

else

    [~,iC] = ...
        min(P3);

    S.phiC = ...
        phi3(iC);

end

%% ------------------------------------------------------------------------
% Final resistance optimum
% -------------------------------------------------------------------------

if S.resistanceFlat

    S.phiR = ...
        NaN;

else

    [~,iR] = ...
        min(R3);

    S.phiR = ...
        phi3(iR);

end

%% ------------------------------------------------------------------------
% Separation and penalties
% -------------------------------------------------------------------------

if ~S.controlFlat && ...
   ~S.resistanceFlat

    S.dCR = ...
        circular_distance( ...
        S.phiC,S.phiR,p.T);

    [PC,RC] = ...
        evaluate_phase( ...
        S.phiC,z0,p);

    [PR,RR] = ...
        evaluate_phase( ...
        S.phiR,z0,p);

    S.controlPenalty = ...
        100*(PR-PC)/PC;

    S.resistancePenaltyPP = ...
        100*(RC-RR);

else

    S.dCR = ...
        NaN;

    S.controlPenalty = ...
        NaN;

    S.resistancePenaltyPP = ...
        NaN;

end

end


%% ========================================================================
% PERIODIC UNTREATED RESIDENT
% ========================================================================

function [z,nYears,lastError] = periodic_resident(p)

Nv = ...
    p.m*p.Nh;

Ih = ...
    0.01*p.Nh;

Eh = ...
    0.005*p.Nh;

Sh = ...
    p.Nh-Eh-Ih;

Iv = ...
    0.01*Nv;

Ev = ...
    0.005*Nv;

Sv = ...
    Nv-Ev-Iv;

z = ...
    [Sh;Eh;Ih;Sv;Ev;Iv];

lastError = ...
    Inf;

for nYears = 1:p.maxBurnYears

    [~,Z] = ode45( ...
        @(t,x) untreated_rhs(t,x,p), ...
        [0 p.T], ...
        z, ...
        p.opts6);

    zNew = ...
        Z(end,:)';

    lastError = ...
        norm(zNew-z,inf) / ...
        max(1,norm(z,inf));

    z = ...
        zNew;

    if lastError < p.periodicTol

        break;

    end

end

end


%% ========================================================================
% UNTREATED RESIDENT AT CALENDAR PHASE
% ========================================================================

function zPhi = resident_at_phase(phi,z0,p)

if phi == 0

    zPhi = ...
        z0;

    return;

end

[~,Z] = ode45( ...
    @(t,z) untreated_rhs(t,z,p), ...
    [0 phi], ...
    z0, ...
    p.opts6);

zPhi = ...
    Z(end,:)';

end


%% ========================================================================
% EVALUATE ONE PHASE
% ========================================================================

function [Pbar,Rpool] = evaluate_phase(phi,z0,p)

zPhi = ...
    resident_at_phase( ...
    phi,z0,p);

y0 = ...
    seed_resistance( ...
    zPhi,p);

tStart = ...
    phi;

tEnd = ...
    phi+p.years*p.T;

[t,Y] = ...
    simulate_smc( ...
    y0,p,tStart,tEnd,phi);

mask = ...
    t >= ...
    tEnd-p.evalYears*p.T;

tt = ...
    t(mask);

IhS = ...
    Y(mask,3);

IhR = ...
    Y(mask,5);

%% Malaria prevalence

Prev = ...
    (IhS+IhR)/p.Nh;

Pbar = ...
    trapz(tt,Prev) / ...
    (tt(end)-tt(1));

%% Pooled resistance

den = ...
    trapz( ...
    tt, ...
    IhS+IhR);

if den > 1e-12

    Rpool = ...
        trapz(tt,IhR)/den;

else

    Rpool = ...
        0;

end

end


%% ========================================================================
% RESISTANCE SEEDING
% ========================================================================

function y = seed_resistance(z,p)

y = ...
    zeros(11,1);

y(1) = z(1);

y(2) = ...
    (1-p.pR0)*z(2);

y(3) = ...
    (1-p.pR0)*z(3);

y(4) = ...
    p.pR0*z(2);

y(5) = ...
    p.pR0*z(3);

y(6) = 0;

y(7) = z(4);

y(8) = ...
    (1-p.pR0)*z(5);

y(9) = ...
    (1-p.pR0)*z(6);

y(10) = ...
    p.pR0*z(5);

y(11) = ...
    p.pR0*z(6);

end


%% ========================================================================
% UNTREATED MODEL
% ========================================================================

function dz = untreated_rhs(t,z,p)

Sh = z(1);
Eh = z(2);
Ih = z(3);

Sv = z(4);
Ev = z(5);
Iv = z(6);

lambda_h = ...
    p.a*p.bS*Iv/p.Nh;

lambda_v = ...
    p.a*p.cS*Ih/p.Nh;

Lambda_v = ...
    p.LambdaBar * ...
    (1 + ...
    p.xi_v*cos( ...
    2*pi*(t-p.phi_v)/p.T));

Lambda_v = ...
    max(Lambda_v,0);

dSh = ...
    p.mu_h*p.Nh ...
    +p.gammaS*Ih ...
    -lambda_h*Sh ...
    -p.mu_h*Sh;

dEh = ...
    lambda_h*Sh ...
    -(p.sigma_h+p.mu_h)*Eh;

dIh = ...
    p.sigma_h*Eh ...
    -(p.gammaS+p.mu_h)*Ih;

dSv = ...
    Lambda_v ...
    -lambda_v*Sv ...
    -p.mu_v*Sv;

dEv = ...
    lambda_v*Sv ...
    -(p.sigma_v+p.mu_v)*Ev;

dIv = ...
    p.sigma_v*Ev ...
    -p.mu_v*Iv;

dz = ...
    [dSh;dEh;dIh; ...
     dSv;dEv;dIv];

end


%% ========================================================================
% FULL MODEL
% ========================================================================

function dy = full_rhs(t,y,p)

Sh = y(1);

EhS = y(2);
IhS = y(3);

EhR = y(4);
IhR = y(5);

Ph = y(6);

Sv = y(7);

EvS = y(8);
IvS = y(9);

EvR = y(10);
IvR = y(11);

%% Forces of infection

lhS = ...
    p.a*p.bS*IvS/p.Nh;

lhR = ...
    p.a*p.bR*IvR/p.Nh;

lvS = ...
    p.a*p.cS*IhS/p.Nh;

lvR = ...
    p.a*p.cR*IhR/p.Nh;

%% Mosquito recruitment

Lv = ...
    p.LambdaBar * ...
    (1 + ...
    p.xi_v*cos( ...
    2*pi*(t-p.phi_v)/p.T));

Lv = ...
    max(Lv,0);

%% Humans

dSh = ...
    p.mu_h*p.Nh ...
    +p.gammaS*IhS ...
    +p.gammaR*IhR ...
    +p.omega*Ph ...
    -lhS*Sh ...
    -lhR*Sh ...
    -p.mu_h*Sh;

dEhS = ...
    lhS*Sh ...
    +(1-p.epsS)*lhS*Ph ...
    -(p.sigma_h+p.mu_h)*EhS;

dIhS = ...
    p.sigma_h*EhS ...
    -(p.gammaS+p.mu_h)*IhS;

dEhR = ...
    lhR*Sh ...
    +(1-p.epsR)*lhR*Ph ...
    -(p.sigma_h+p.mu_h)*EhR;

dIhR = ...
    p.sigma_h*EhR ...
    -(p.gammaR+p.mu_h)*IhR;

dPh = ...
    -(1-p.epsS)*lhS*Ph ...
    -(1-p.epsR)*lhR*Ph ...
    -(p.omega+p.mu_h)*Ph;

%% Mosquitoes

dSv = ...
    Lv ...
    -lvS*Sv ...
    -lvR*Sv ...
    -p.mu_v*Sv;

dEvS = ...
    lvS*Sv ...
    -(p.sigma_v+p.mu_v)*EvS;

dIvS = ...
    p.sigma_v*EvS ...
    -p.mu_v*IvS;

dEvR = ...
    lvR*Sv ...
    -(p.sigma_v+p.mu_v)*EvR;

dIvR = ...
    p.sigma_v*EvR ...
    -p.mu_v*IvR;

dy = ...
    [dSh;dEhS;dIhS; ...
     dEhR;dIhR;dPh; ...
     dSv;dEvS;dIvS; ...
     dEvR;dIvR];

end


%% ========================================================================
% SMC PULSE
% ========================================================================

function y = apply_pulse(y,p)

Sh = ...
    y(1);

IhS = ...
    y(3);

IhR = ...
    y(5);

Ph = ...
    y(6);

treated = ...
    p.q*Sh;

clearS = ...
    p.q*p.eS*IhS;

clearR = ...
    p.q*p.eR*IhR;

y(1) = ...
    Sh-treated;

y(3) = ...
    IhS-clearS;

y(5) = ...
    IhR-clearR;

y(6) = ...
    Ph ...
    +treated ...
    +clearS ...
    +clearR;

end


%% ========================================================================
% FAST SMC SIMULATION
% ========================================================================

function [tout,Yout] = ...
    simulate_smc(y0,p,tStart,tEnd,phi)

pulse = ...
    zeros(p.years*p.K,1);

idx = 0;

for n = 0:(p.years-1)

    for k = 0:(p.K-1)

        idx = ...
            idx+1;

        pulse(idx) = ...
            phi ...
            +n*p.T ...
            +k*p.tau;

    end

end

tcur = ...
    tStart;

ycur = ...
    y0;

tout = ...
    tStart;

Yout = ...
    y0';

for j = 1:length(pulse)

    tp = ...
        pulse(j);

    if tp > tcur+1e-12

        [ts,ys] = ode45( ...
            @(t,y) full_rhs(t,y,p), ...
            [tcur tp], ...
            ycur, ...
            p.opts11);

        tout = ...
            [tout;ts(2:end)]; %#ok<AGROW>

        Yout = ...
            [Yout;ys(2:end,:)]; %#ok<AGROW>

        ycur = ...
            ys(end,:)';

    end

    ycur = ...
        apply_pulse(ycur,p);

    tout = ...
        [tout;tp]; %#ok<AGROW>

    Yout = ...
        [Yout;ycur']; %#ok<AGROW>

    tcur = ...
        tp;

end

if tcur < tEnd

    [ts,ys] = ode45( ...
        @(t,y) full_rhs(t,y,p), ...
        [tcur tEnd], ...
        ycur, ...
        p.opts11);

    tout = ...
        [tout;ts(2:end)];

    Yout = ...
        [Yout;ys(2:end,:)];

end

end


%% ========================================================================
% CIRCULAR DISTANCE
% ========================================================================

function d = circular_distance(a,b,T)

x = ...
    abs(a-b);

d = ...
    min(x,T-x);

end


%% ========================================================================
% CLEAN PLOTTING FUNCTION
%
% NaNs are left as gaps.
% No NI annotations are displayed.
% ========================================================================

function plot_optima(ax,x,R,xLabel,panelTitle)

hold(ax,'on');

%% Malaria-control optimum

plot( ...
    ax, ...
    x, ...
    R.phiC, ...
    'o-', ...
    'LineWidth',1.6, ...
    'MarkerSize',6, ...
    'MarkerFaceColor','w');

%% Resistance optimum
%
% MATLAB automatically leaves gaps where phiR = NaN.

plot( ...
    ax, ...
    x, ...
    R.phiR, ...
    's-', ...
    'LineWidth',1.6, ...
    'MarkerSize',6, ...
    'MarkerFaceColor','w');

xlabel( ...
    ax, ...
    xLabel);

ylabel( ...
    ax, ...
    'Optimal starting phase (day)');

title( ...
    ax, ...
    panelTitle, ...
    'FontWeight','bold');

hold(ax,'off');

end