%% ========================================================================
% CODE 3 -- FINAL STANDALONE CENTRAL PHASE OPTIMIZATION
%
% Computes:
%   phi_C* : phase minimizing long-term malaria burden
%   phi_R* : phase minimizing established resistant frequency
%   d_CR   : circular separation
%   cross-optimum penalties and Pareto phases
%
% FINAL REFERENCE
%   gamma = 1/20 day^-1, xi_v = 0.40
%   fR = 0.95, epsR = 0.50
%   q = 0.70, K = 4, tau = 30 d
% ========================================================================

clear; clc; close all;

%% Parameters
p.Nh=1e5; p.a=0.33; p.bS=0.20; p.cS=0.40; p.m=2.0;
p.sigma_h=1/12; p.sigma_v=1/10;
p.gammaS=1/20; p.gammaR=1/20;
p.mu_h=1/(65*365); p.mu_v=0.0714;

p.T=365; p.xi_v=0.40; p.phi_v=0;
p.LambdaBar=p.m*p.Nh*p.mu_v;

p.fR=0.95;
p.bR=p.fR*p.bS;
p.cR=p.fR*p.cS;

p.K=4; p.tau=30; p.q=0.70;
p.eS=0.98; p.eR=0.90;
p.epsS=0.90; p.epsR=0.50;
p.DP=35; p.omega=1/p.DP;

p.pR0=0.01;
p.burnYears=80; p.years=20; p.evalYears=5;
p.resistanceFlatTolerancePP=0.10;

p.opts6=odeset('RelTol',1e-9,'AbsTol',1e-11,...
               'MaxStep',1,'NonNegative',1:6);
p.opts11=odeset('RelTol',1e-9,'AbsTol',1e-11,...
                'MaxStep',1,'NonNegative',1:11);

R0S=calc_R0(p,p.bS,p.cS,p.gammaS);
R0R=calc_R0(p,p.bR,p.cR,p.gammaR);

fprintf('\n============================================================\n');
fprintf('CODE 3: CENTRAL SEASONAL SMC PHASE OPTIMIZATION\n');
fprintf('============================================================\n');
fprintf('R0S = %.6f\n',R0S);
fprintf('R0R = %.6f\n',R0R);
fprintf('Relative resistant fitness fR = %.2f\n',p.fR);
fprintf('Seasonality xi_v = %.2f\n',p.xi_v);
fprintf('SMC q = %.2f, K = %d, tau = %.0f d\n',p.q,p.K,p.tau);

%% Untreated periodic resident at phase zero
fprintf('\nRunning untreated seasonal burn-in...\n');
zPhase0=get_phase_zero_resident(p);
P0=evaluate_untreated_reference(zPhase0,p);
fprintf('Untreated reference mean prevalence = %.3f %%\n',100*P0);

%% Full daily phase sweep
phiAll=0:(p.T-1);
nPhi=numel(phiAll);
Pall=nan(1,nPhi); Rall=nan(1,nPhi);
TerminalR=nan(1,nPhi); pulseCounts=zeros(1,nPhi);
humanErrors=nan(1,nPhi); seedErrors=nan(1,nPhi);

fprintf('\nRunning complete one-day phase sweep...\n');

for i=1:nPhi
    phi=phiAll(i);
    zPhi=resident_at_phase(phi,zPhase0,p);
    ySeed=seed_resistance(zPhi,p);
    initialRes=ySeed(5)/(ySeed(3)+ySeed(5));
    seedErrors(i)=abs(initialRes-p.pR0);

    [Pall(i),Rall(i),TerminalR(i),pulseCounts(i),humanErrors(i)] = ...
        evaluate_phase(phi,ySeed,p);

    if mod(phi,20)==0 || phi==364
        fprintf('Phase %3d | burden = %6.3f%% | pooled resistance = %6.3f%% | terminal = %6.3f%%\n',...
            phi,100*Pall(i),100*Rall(i),100*TerminalR(i));
    end
end

expectedPulses=p.years*p.K;
if any(pulseCounts~=expectedPulses)
    error('Pulse-count error: every phase must receive %d pulses.',expectedPulses);
end

%% Optima
[Pmin,iC]=min(Pall); phiC=phiAll(iC);
Pmax=max(Pall);

Rmin=min(Rall); Rmax=max(Rall);
resistanceRangePP=100*(Rmax-Rmin);
resistanceMeaningful=resistanceRangePP>=p.resistanceFlatTolerancePP;

if resistanceMeaningful
    [Rmin,iR]=min(Rall); phiR=phiAll(iR);
    dCR=circ_dist(phiC,phiR,p.T);
else
    iR=[]; phiR=NaN; dCR=NaN;
end

EC_phiC=100*(1-Pmin/P0);

if resistanceMeaningful
    Pi_C_from_R=100*(Pall(iR)-Pmin)/Pmin;
    Pi_R_from_C_pp=100*(Rall(iC)-Rmin);
    Pi_R_from_C_rel=100*(Rall(iC)-Rmin)/Rmin;
else
    Pi_C_from_R=NaN; Pi_R_from_C_pp=NaN; Pi_R_from_C_rel=NaN;
end

%% Pareto
ND=false(size(phiAll));
if resistanceMeaningful
    ND(:)=true;
    for i=1:nPhi
        dominated=any(Pall<=Pall(i) & Rall<=Rall(i) & ...
                     (Pall<Pall(i) | Rall<Rall(i)));
        if dominated, ND(i)=false; end
    end
end

%% Output
fprintf('\n============================================================\n');
fprintf('CENTRAL PHASE-OPTIMIZATION RESULTS\n');
fprintf('============================================================\n');
fprintf('Control optimum phi_C* = %d days\n',phiC);
fprintf('Malaria burden range = %.3f -- %.3f %%\n',100*min(Pall),100*max(Pall));
fprintf('Malaria phase range = %.3f pp\n',100*(max(Pall)-min(Pall)));
fprintf('SMC effect at phi_C* = %.3f %%\n',EC_phiC);
fprintf('Established resistance range = %.3f -- %.3f %%\n',100*Rmin,100*Rmax);
fprintf('Resistance phase range Delta_R = %.3f pp\n',resistanceRangePP);

if resistanceMeaningful
    fprintf('Resistance optimum phi_R* = %d days\n',phiR);
    fprintf('Circular separation d_CR = %.1f days\n',dCR);
    fprintf('Control penalty at phi_R* = %.4f %%\n',Pi_C_from_R);
    fprintf('Resistance penalty at phi_C* = %.4f pp\n',Pi_R_from_C_pp);
    fprintf('Relative resistance penalty at phi_C* = %.4f %%\n',Pi_R_from_C_rel);
    fprintf('Number of Pareto phases = %d\n',sum(ND));
else
    fprintf('Resistance objective effectively flat: no meaningful phi_R* declared.\n');
end

fprintf('Expected pulses per phase = %d\n',expectedPulses);
fprintf('Maximum human conservation error = %.4e\n',max(humanErrors));
fprintf('Maximum resistance-seeding error = %.4e\n',max(seedErrors));
fprintf('============================================================\n');

PhaseResults=table(phiAll',100*Pall',100*Rall',100*TerminalR',...
                   pulseCounts',humanErrors',...
    'VariableNames',{'SMC_start_day','MeanMalariaPrev_percent',...
    'PooledResistance_percent','TerminalResistance_percent',...
    'PulseCount','HumanConservationError'});

%% Figure
figure('Color','w','Position',[60 100 1450 480]);
TL=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

nexttile;
plot(phiAll,100*Pall,'LineWidth',1.7); hold on;
plot(phiC,100*Pall(iC),'kp','MarkerSize',12,'MarkerFaceColor','k');
xline(phiC,'--');
xlabel('SMC starting phase, \phi (day)');
ylabel('Mean malaria prevalence (%)');
title('(a) Long-term malaria burden');
xlim([0 364]); grid on; box on;
text(phiC+5,100*Pall(iC),sprintf('\\phi_C^*=%d d',phiC),...
    'VerticalAlignment','bottom','FontSize',9); hold off;

nexttile;
plot(phiAll,100*Rall,'LineWidth',1.7); hold on;
if resistanceMeaningful
    plot(phiR,100*Rall(iR),'kp','MarkerSize',12,'MarkerFaceColor','k');
    xline(phiR,'--');
    text(phiR+5,100*Rall(iR),sprintf('\\phi_R^*=%d d',phiR),...
        'VerticalAlignment','bottom','FontSize',9);
end
xlabel('SMC starting phase, \phi (day)');
ylabel('Pooled resistant frequency (%)');
title('(b) Established resistance');
xlim([0 364]); grid on; box on; hold off;

nexttile;
plot(100*Pall,100*Rall,'o','MarkerSize',4); hold on;
plot(100*Pall(iC),100*Rall(iC),'kp','MarkerSize',12,'MarkerFaceColor','k');
if resistanceMeaningful
    plot(100*Pall(iR),100*Rall(iR),'ks','MarkerSize',9,'MarkerFaceColor','k');
    plot(100*Pall(ND),100*Rall(ND),'o','MarkerSize',7,'LineWidth',1.2);
    legend('All phases','\phi_C^*','\phi_R^*','Pareto phases','Location','best');
else
    legend('All phases','\phi_C^*','Location','best');
end
xlabel('Mean malaria prevalence (%)');
ylabel('Pooled resistant frequency (%)');
title('(c) Control--resistance trade-off');
grid on; box on; hold off;

if resistanceMeaningful
    %title(TL,sprintf('Seasonal SMC timing: \\phi_C^*=%d d, \\phi_R^*=%d d, d_{CR}=%.0f d',...
       % phiC,phiR,dCR),'FontSize',13,'FontWeight','bold');
else
    title(TL,sprintf('Seasonal SMC timing: \\phi_C^*=%d d; resistance phase-insensitive',...
        phiC),'FontSize',13,'FontWeight','bold');
end

%% Local functions
function R0=calc_R0(p,b,c,gamma)
R0=sqrt((p.m*p.a^2*b*c/((gamma+p.mu_h)*p.mu_v))*...
        (p.sigma_h/(p.sigma_h+p.mu_h))*...
        (p.sigma_v/(p.sigma_v+p.mu_v)));
end

function d=circ_dist(a,b,T)
x=abs(a-b); d=min(x,T-x);
end

function z0=get_phase_zero_resident(p)
Nv=p.m*p.Nh;
Ih=.01*p.Nh; Eh=.005*p.Nh; Sh=p.Nh-Eh-Ih;
Iv=.01*Nv; Ev=.005*Nv; Sv=Nv-Ev-Iv;
zinit=[Sh;Eh;Ih;Sv;Ev;Iv];
[~,Z]=ode15s(@(t,z) untreated_rhs(t,z,p),...
             [0 p.burnYears*p.T],zinit,p.opts6);
z0=Z(end,:)';
end

function zPhi=resident_at_phase(phi,z0,p)
if phi==0, zPhi=z0; return; end
[~,Z]=ode15s(@(t,z) untreated_rhs(t,z,p),[0 phi],z0,p.opts6);
zPhi=Z(end,:)';
end

function y=seed_resistance(z,p)
y=zeros(11,1);
y(1)=z(1);
y(2)=(1-p.pR0)*z(2); y(4)=p.pR0*z(2);
y(3)=(1-p.pR0)*z(3); y(5)=p.pR0*z(3);
y(6)=0;
y(7)=z(4);
y(8)=(1-p.pR0)*z(5); y(10)=p.pR0*z(5);
y(9)=(1-p.pR0)*z(6); y(11)=p.pR0*z(6);
end

function P0=evaluate_untreated_reference(z0,p)
tEnd=5*p.T;
[t,Z]=ode15s(@(t,z) untreated_rhs(t,z,p),[0 tEnd],z0,p.opts6);
P=Z(:,3)/p.Nh;
P0=trapz(t,P)/(t(end)-t(1));
end

function [Pbar,Rpool,Rterminal,pulseCount,humanError]=evaluate_phase(phi,y0,p)
tStart=phi; tEnd=phi+p.years*p.T;
[t,Y,pulseCount]=simulate_smc(y0,p,tStart,tEnd,phi);
mask=t>=tEnd-p.evalYears*p.T; tt=t(mask);
IhS=Y(mask,3); IhR=Y(mask,5);
P=(IhS+IhR)/p.Nh;
Pbar=trapz(tt,P)/(tt(end)-tt(1));
den=trapz(tt,IhS+IhR);
if den>1e-12, Rpool=trapz(tt,IhR)/den; else, Rpool=NaN; end
d=Y(end,3)+Y(end,5);
if d>1e-12, Rterminal=Y(end,5)/d; else, Rterminal=NaN; end
humanError=max(abs(sum(Y(:,1:6),2)-p.Nh));
end

function dz=untreated_rhs(t,z,p)
Sh=z(1); Eh=z(2); Ih=z(3); Sv=z(4); Ev=z(5); Iv=z(6);
lh=p.a*p.bS*Iv/p.Nh; lv=p.a*p.cS*Ih/p.Nh;
Lv=p.LambdaBar*(1+p.xi_v*cos(2*pi*(t-p.phi_v)/p.T)); Lv=max(Lv,0);
dz=[p.mu_h*p.Nh+p.gammaS*Ih-lh*Sh-p.mu_h*Sh;
    lh*Sh-(p.sigma_h+p.mu_h)*Eh;
    p.sigma_h*Eh-(p.gammaS+p.mu_h)*Ih;
    Lv-lv*Sv-p.mu_v*Sv;
    lv*Sv-(p.sigma_v+p.mu_v)*Ev;
    p.sigma_v*Ev-p.mu_v*Iv];
end

function dy=full_rhs(t,y,p)
Sh=y(1); EhS=y(2); IhS=y(3); EhR=y(4); IhR=y(5); Ph=y(6);
Sv=y(7); EvS=y(8); IvS=y(9); EvR=y(10); IvR=y(11);
lhS=p.a*p.bS*IvS/p.Nh; lhR=p.a*p.bR*IvR/p.Nh;
lvS=p.a*p.cS*IhS/p.Nh; lvR=p.a*p.cR*IhR/p.Nh;
Lv=p.LambdaBar*(1+p.xi_v*cos(2*pi*(t-p.phi_v)/p.T)); Lv=max(Lv,0);

dSh=p.mu_h*p.Nh+p.gammaS*IhS+p.gammaR*IhR+p.omega*Ph-lhS*Sh-lhR*Sh-p.mu_h*Sh;
dEhS=lhS*Sh+(1-p.epsS)*lhS*Ph-(p.sigma_h+p.mu_h)*EhS;
dIhS=p.sigma_h*EhS-(p.gammaS+p.mu_h)*IhS;
dEhR=lhR*Sh+(1-p.epsR)*lhR*Ph-(p.sigma_h+p.mu_h)*EhR;
dIhR=p.sigma_h*EhR-(p.gammaR+p.mu_h)*IhR;
dPh=-(1-p.epsS)*lhS*Ph-(1-p.epsR)*lhR*Ph-(p.omega+p.mu_h)*Ph;
dSv=Lv-lvS*Sv-lvR*Sv-p.mu_v*Sv;
dEvS=lvS*Sv-(p.sigma_v+p.mu_v)*EvS;
dIvS=p.sigma_v*EvS-p.mu_v*IvS;
dEvR=lvR*Sv-(p.sigma_v+p.mu_v)*EvR;
dIvR=p.sigma_v*EvR-p.mu_v*IvR;
dy=[dSh;dEhS;dIhS;dEhR;dIhR;dPh;dSv;dEvS;dIvS;dEvR;dIvR];
end

function y=apply_pulse(y,p)
Sh=y(1); IhS=y(3); IhR=y(5); Ph=y(6);
treated=p.q*Sh; clearS=p.q*p.eS*IhS; clearR=p.q*p.eR*IhR;
y(1)=Sh-treated; y(3)=IhS-clearS; y(5)=IhR-clearR;
y(6)=Ph+treated+clearS+clearR;
end

function [tout,Yout,pulseCount]=simulate_smc(y0,p,tStart,tEnd,phi)
pulse=zeros(p.years*p.K,1); idx=0;
for n=0:(p.years-1)
    for k=0:(p.K-1)
        idx=idx+1; pulse(idx)=phi+n*p.T+k*p.tau;
    end
end
pulseCount=length(pulse);
if any(diff(pulse)<=0), error('SMC pulse sequence not chronological.'); end
if pulse(end)>=tEnd, error('Last SMC pulse outside intervention interval.'); end

tcur=tStart; ycur=y0; tout=tStart; Yout=y0';
for j=1:pulseCount
    tp=pulse(j);
    if tp>tcur+1e-12
        [ts,ys]=ode15s(@(t,y) full_rhs(t,y,p),[tcur tp],ycur,p.opts11);
        tout=[tout;ts(2:end)]; Yout=[Yout;ys(2:end,:)]; %#ok<AGROW>
        ycur=ys(end,:)';
    end
    H0=sum(ycur(1:6));
    ycur=apply_pulse(ycur,p);
    if abs(sum(ycur(1:6))-H0)>1e-6, error('Human population not conserved at pulse.'); end
    tout=[tout;tp]; Yout=[Yout;ycur']; %#ok<AGROW>
    tcur=tp;
end
if tcur<tEnd
    [ts,ys]=ode15s(@(t,y) full_rhs(t,y,p),[tcur tEnd],ycur,p.opts11);
    tout=[tout;ts(2:end)]; Yout=[Yout;ys(2:end,:)];
end
end
