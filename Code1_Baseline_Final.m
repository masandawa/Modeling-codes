%% ========================================================================
% CODE 1 -- FINAL STANDALONE BASELINE
% UNTREATED SEASONAL MALARIA DYNAMICS
%
% Purpose
%   1. Verify the final reference transmission parameterization.
%   2. Compute R0S.
%   3. Obtain the untreated periodic seasonal reference state.
%   4. Report annual malaria prevalence and seasonal peak timing.
%
% FINAL REFERENCE SETTING
%   m        = 2
%   gamma    = 1/20 day^-1
%   xi_v     = 0.40
%   R0S      ~ 1.686
% ========================================================================

clear; clc; close all;

%% Parameters
p.Nh = 1e5;
p.a  = 0.33;
p.bS = 0.20;
p.cS = 0.40;
p.m  = 2.0;

p.sigma_h = 1/12;
p.sigma_v = 1/10;
p.gammaS  = 1/20;

p.mu_h = 1/(65*365);
p.mu_v = 0.0714;

p.T     = 365;
p.xi_v  = 0.40;
p.phi_v = 0;

p.LambdaBar = p.m*p.Nh*p.mu_v;

opts = odeset('RelTol',1e-9,'AbsTol',1e-11,...
              'MaxStep',1.0,'NonNegative',1:6);

%% R0
R0S = sqrt( ...
    (p.m*p.a^2*p.bS*p.cS / ((p.gammaS+p.mu_h)*p.mu_v)) ...
    *(p.sigma_h/(p.sigma_h+p.mu_h)) ...
    *(p.sigma_v/(p.sigma_v+p.mu_v)) );

fprintf('\n============================================================\n');
fprintf('CODE 1: UNTREATED SEASONAL BASELINE\n');
fprintf('============================================================\n');
fprintf('R0S = %.6f\n',R0S);
fprintf('Seasonality xi_v = %.2f\n',p.xi_v);
fprintf('Effective infectious duration = %.1f days\n',1/p.gammaS);

%% Initial conditions
Nv0 = p.m*p.Nh;

Ih0 = 0.01*p.Nh;
Eh0 = 0.005*p.Nh;
Sh0 = p.Nh-Eh0-Ih0;

Iv0 = 0.01*Nv0;
Ev0 = 0.005*Nv0;
Sv0 = Nv0-Ev0-Iv0;

z0 = [Sh0;Eh0;Ih0;Sv0;Ev0;Iv0];

%% Burn-in
burnYears = 80;
fprintf('\nRunning %d-year untreated burn-in...\n',burnYears);

[~,Zburn] = ode15s(@(t,z) rhs(t,z,p),...
                   [0 burnYears*p.T],z0,opts);

zStar = Zburn(end,:)';

%% Final complete year on the periodic attractor
t0 = burnYears*p.T;
t1 = t0+p.T;

tEval = (t0:1:t1)';

[t,Z] = ode15s(@(t,z) rhs(t,z,p),tEval,zStar,opts);

Sh = Z(:,1); Eh = Z(:,2); Ih = Z(:,3);
Sv = Z(:,4); Ev = Z(:,5); Iv = Z(:,6);

NhCheck = Sh+Eh+Ih;
Nv = Sv+Ev+Iv;

humanPrev = Ih/p.Nh;
mosqPrev  = Iv./Nv;
mosqRatio = Nv/p.Nh;

day = t-t0;

meanPrev = trapz(t,humanPrev)/(t(end)-t(1));
minPrev  = min(humanPrev);
maxPrev  = max(humanPrev);
ampPP    = 100*(maxPrev-minPrev);

[~,iNv] = max(Nv);
[~,iIv] = max(mosqPrev);
[~,iIh] = max(humanPrev);

dayNv = day(iNv);
dayIv = day(iIv);
dayIh = day(iIh);

lagNvToHuman = mod(dayIh-dayNv,p.T);

humanError = max(abs(NhCheck-p.Nh));

fprintf('\n---------------- BASELINE OUTPUT ----------------\n');
fprintf('Mean annual malaria prevalence = %.3f %%\n',100*meanPrev);
fprintf('Minimum malaria prevalence     = %.3f %%\n',100*minPrev);
fprintf('Maximum malaria prevalence     = %.3f %%\n',100*maxPrev);
fprintf('Seasonal prevalence range      = %.3f pp\n',ampPP);
fprintf('\nMosquito abundance peak day    = %.0f\n',dayNv);
fprintf('Infectious mosquito peak day   = %.0f\n',dayIv);
fprintf('Human malaria peak day         = %.0f\n',dayIh);
fprintf('Mosquito abundance -> human lag= %.0f days\n',lagNvToHuman);
fprintf('Maximum human population error = %.4e\n',humanError);
fprintf('-------------------------------------------------\n');

%% Figure
figure('Color','w','Position',[100 70 1000 780]);

subplot(3,1,1)
plot(day,mosqRatio,'LineWidth',1.6);
xlabel('Day of year'); ylabel('N_v/N_h');
title('(a) Seasonal mosquito abundance');
xlim([0 365]); grid on; box on;

subplot(3,1,2)
plot(day,100*mosqPrev,'LineWidth',1.6);
xlabel('Day of year'); ylabel('Infectious mosquitoes (%)');
title('(b) Infectious-mosquito prevalence');
xlim([0 365]); grid on; box on;

subplot(3,1,3)
plot(day,100*humanPrev,'LineWidth',1.6);
xlabel('Day of year'); ylabel('Human malaria prevalence (%)');
title('(c) Human infectious prevalence');
xlim([0 365]); grid on; box on;

%sgtitle(sprintf('Untreated seasonal reference: R_0 = %.2f, \\xi_v = %.2f',R0S,p.xi_v));

%% Local function
function dz = rhs(t,z,p)
Sh=z(1); Eh=z(2); Ih=z(3);
Sv=z(4); Ev=z(5); Iv=z(6);

lambda_h = p.a*p.bS*Iv/p.Nh;
lambda_v = p.a*p.cS*Ih/p.Nh;

Lambda_v = p.LambdaBar*(1+p.xi_v*cos(2*pi*(t-p.phi_v)/p.T));
Lambda_v = max(Lambda_v,0);

dSh = p.mu_h*p.Nh+p.gammaS*Ih-lambda_h*Sh-p.mu_h*Sh;
dEh = lambda_h*Sh-(p.sigma_h+p.mu_h)*Eh;
dIh = p.sigma_h*Eh-(p.gammaS+p.mu_h)*Ih;

dSv = Lambda_v-lambda_v*Sv-p.mu_v*Sv;
dEv = lambda_v*Sv-(p.sigma_v+p.mu_v)*Ev;
dIv = p.sigma_v*Ev-p.mu_v*Iv;

dz=[dSh;dEh;dIh;dSv;dEv;dIv];
end
