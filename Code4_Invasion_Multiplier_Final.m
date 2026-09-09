%% ========================================================================
% CODE 4 -- FINAL STANDALONE
% PERIODIC RARE-RESISTANCE INVASION MULTIPLIER
%
% Full daily phase sweep (0:364) using Floquet/monodromy analysis.
%
% FINAL REFERENCE
%   R0S ~ 1.686, R0R ~ 1.602
%   fR = 0.95, xi_v = 0.40
%   q = 0.70, K = 4, tau = 30 d
%
% Code 3 reference optima:
%   phi_C* = 102 d
%   phi_R* = 118 d
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

p.periodicTol=1e-9;
p.maxPeriodicYears=200;

p.optsSensitive=odeset('RelTol',1e-9,'AbsTol',1e-11,...
                        'MaxStep',0.5,'NonNegative',1:7);
p.optsAugmented=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',0.5);

R0S=calc_R0(p,p.bS,p.cS,p.gammaS);
R0R=calc_R0(p,p.bR,p.cR,p.gammaR);

fprintf('\n============================================================\n');
fprintf('CODE 4: PERIODIC RARE-RESISTANCE INVASION\n');
fprintf('============================================================\n');
fprintf('R0S = %.6f\n',R0S);
fprintf('R0R = %.6f\n',R0R);
fprintf('Relative resistant fitness fR = %.2f\n',p.fR);
fprintf('Seasonality xi_v = %.2f\n',p.xi_v);
fprintf('SMC q = %.2f, K = %d, tau = %.0f d\n',p.q,p.K,p.tau);

yInitial=initial_sensitive_state(p);

%% Full daily phase sweep
phiAll=0:(p.T-1);
Mall=nan(size(phiAll));
ConvYears=nan(size(phiAll));
PeriodicErr=nan(size(phiAll));

fprintf('\nRunning complete one-day invasion sweep...\n');

for i=1:length(phiAll)
    phi=phiAll(i);
    [Mall(i),ConvYears(i),PeriodicErr(i)] = invasion_multiplier(phi,yInitial,p);

    if mod(phi,20)==0 || phi==364
        fprintf('Phase %3d | M_R = %.6f | periodic years = %d | error = %.3e\n',...
            phi,Mall(i),ConvYears(i),PeriodicErr(i));
    end
end

[Mmin,iI]=min(Mall); Mmax=max(Mall); phiI=phiAll(iI);
Mrange=Mmax-Mmin;
relativeMrange=100*Mrange/Mmin;

phiC_ref=102; phiR_ref=118;
dCI=circ_dist(phiC_ref,phiI,p.T);
dRI=circ_dist(phiR_ref,phiI,p.T);

fprintf('\n============================================================\n');
fprintf('RARE-RESISTANCE INVASION RESULTS\n');
fprintf('============================================================\n');
fprintf('Invasion optimum phi_I* = %d days\n',phiI);
fprintf('Minimum M_R = %.6f\n',Mmin);
fprintf('Maximum M_R = %.6f\n',Mmax);
fprintf('Absolute M_R phase range = %.6f\n',Mrange);
fprintf('Relative M_R phase range = %.3f %%\n',relativeMrange);

if all(Mall>1)
    fprintf('Rare resistance can invade at EVERY phase.\n');
elseif all(Mall<1)
    fprintf('Rare resistance declines at EVERY phase.\n');
else
    fprintf('M_R crosses 1: invasion depends on SMC timing.\n');
end

fprintf('phi_C* = %d d, phi_R* = %d d, phi_I* = %d d\n',...
    phiC_ref,phiR_ref,phiI);
fprintf('d_CI = %.1f d, d_RI = %.1f d\n',dCI,dRI);
fprintf('Maximum periodic-orbit error = %.3e\n',max(PeriodicErr));
fprintf('Maximum years to convergence = %.0f\n',max(ConvYears));
fprintf('============================================================\n');

InvasionResults=table(phiAll',Mall',ConvYears',PeriodicErr',...
    'VariableNames',{'SMC_start_day','InvasionMultiplier_MR',...
    'PeriodicConvergenceYears','PeriodicConvergenceError'});

%% Clean single-panel figure
figure('Color','w','Position',[100 100 1050 620]);

hM=plot(phiAll,Mall,'-','LineWidth',1.7); hold on;
hThr=yline(1,'--','LineWidth',1.3);
hC=xline(phiC_ref,'--','LineWidth',1.1);
hR=xline(phiR_ref,':','LineWidth',1.4);
hI=plot(phiI,Mmin,'kp','MarkerSize',13,'MarkerFaceColor','k');

text(phiI+7,Mmin+0.012,sprintf('\\phi_I^* = %d d',phiI),...
    'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',10);
text(360,1.008,'M_R = 1','HorizontalAlignment','right',...
    'VerticalAlignment','bottom','FontSize',10);
text(185,1.055,'M_R > 1 for all \phi',...
    'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10);

xlabel('SMC starting phase, \phi (day)');
ylabel('Annual rare-resistance invasion multiplier, M_R');
xlim([0 364]); ylim([0.95 max(Mall)*1.025]);
grid on; box on;
set(gca,'FontSize',11,'LineWidth',0.8);

legend([hM hThr hC hR hI],...
    {'M_R(\phi)','Invasion threshold',...
     '\phi_C^* = 102 d','\phi_R^* = 118 d',...
     sprintf('\\phi_I^* = %d d',phiI)},...
    'Location','northwest','FontSize',9);

%title('Rare-resistance invasion potential across SMC starting phase',...
      %'FontWeight','bold','FontSize',13);
hold off;

%% Local functions
function R0=calc_R0(p,b,c,gamma)
R0=sqrt((p.m*p.a^2*b*c/((gamma+p.mu_h)*p.mu_v))*...
        (p.sigma_h/(p.sigma_h+p.mu_h))*...
        (p.sigma_v/(p.sigma_v+p.mu_v)));
end

function d=circ_dist(a,b,T)
x=abs(a-b); d=min(x,T-x);
end

function y0=initial_sensitive_state(p)
Nv=p.m*p.Nh;
Ih=.01*p.Nh; Eh=.005*p.Nh; Sh=p.Nh-Eh-Ih;
Iv=.01*Nv; Ev=.005*Nv; Sv=Nv-Ev-Iv;
y0=[Sh;Eh;Ih;0;Sv;Ev;Iv];
end

function [M,nYears,lastErr]=invasion_multiplier(phi,yInitial,p)
y=yInitial; lastErr=Inf;
for nYears=1:p.maxPeriodicYears
    yn=sensitive_annual_map(y,phi,p);
    lastErr=norm(yn-y,inf)/max(1,norm(y,inf));
    y=yn;
    if lastErr<p.periodicTol, break; end
end
if nYears==p.maxPeriodicYears && lastErr>=p.periodicTol
    warning('Periodic sensitive resident did not converge for phi=%d.',phi);
end

Phi0=eye(4);
aug0=[y;Phi0(:)];
augT=augmented_annual_map(aug0,phi,p);
PhiT=reshape(augT(8:end),4,4);
M=max(abs(eig(PhiT)));
end

function yT=sensitive_annual_map(y0,phi,p)
pulse=(0:p.K-1)*p.tau;
y=y0; sCur=0;
for j=1:length(pulse)
    sp=pulse(j);
    if sp>sCur+1e-12
        [~,Y]=ode15s(@(s,x) sensitive_rhs(s,x,phi,p),...
                     [sCur sp],y,p.optsSensitive);
        y=Y(end,:)';
    end
    y=sensitive_pulse(y,p);
    sCur=sp;
end
if sCur<p.T
    [~,Y]=ode15s(@(s,x) sensitive_rhs(s,x,phi,p),...
                 [sCur p.T],y,p.optsSensitive);
    y=Y(end,:)';
end
yT=y;
end

function augT=augmented_annual_map(aug0,phi,p)
pulse=(0:p.K-1)*p.tau;
aug=aug0; sCur=0;
for j=1:length(pulse)
    sp=pulse(j);
    if sp>sCur+1e-12
        [~,A]=ode15s(@(s,x) augmented_rhs(s,x,phi,p),...
                     [sCur sp],aug,p.optsAugmented);
        aug=A(end,:)';
    end
    aug=augmented_pulse(aug,p);
    sCur=sp;
end
if sCur<p.T
    [~,A]=ode15s(@(s,x) augmented_rhs(s,x,phi,p),...
                 [sCur p.T],aug,p.optsAugmented);
    aug=A(end,:)';
end
augT=aug;
end

function dy=sensitive_rhs(s,y,phi,p)
tAbs=phi+s;
Sh=y(1); Eh=y(2); Ih=y(3); Ph=y(4);
Sv=y(5); Ev=y(6); Iv=y(7);

lh=p.a*p.bS*Iv/p.Nh;
lv=p.a*p.cS*Ih/p.Nh;
Lv=p.LambdaBar*(1+p.xi_v*cos(2*pi*(tAbs-p.phi_v)/p.T));
Lv=max(Lv,0);

dSh=p.mu_h*p.Nh+p.gammaS*Ih+p.omega*Ph-lh*Sh-p.mu_h*Sh;
dEh=lh*Sh+(1-p.epsS)*lh*Ph-(p.sigma_h+p.mu_h)*Eh;
dIh=p.sigma_h*Eh-(p.gammaS+p.mu_h)*Ih;
dPh=-(1-p.epsS)*lh*Ph-(p.omega+p.mu_h)*Ph;
dSv=Lv-lv*Sv-p.mu_v*Sv;
dEv=lv*Sv-(p.sigma_v+p.mu_v)*Ev;
dIv=p.sigma_v*Ev-p.mu_v*Iv;

dy=[dSh;dEh;dIh;dPh;dSv;dEv;dIv];
end

function y=sensitive_pulse(y,p)
Sh=y(1); Ih=y(3); Ph=y(4);
treated=p.q*Sh;
cleared=p.q*p.eS*Ih;
y(1)=Sh-treated;
y(3)=Ih-cleared;
y(4)=Ph+treated+cleared;
end

function da=augmented_rhs(s,aug,phi,p)
y=aug(1:7);
Phi=reshape(aug(8:end),4,4);
dy=sensitive_rhs(s,y,phi,p);

Sh=y(1); Ph=y(4); Sv=y(5);
kh=p.sigma_h+p.mu_h;
rh=p.gammaR+p.mu_h;
kv=p.sigma_v+p.mu_v;

A=zeros(4);
A(1,1)=-kh;
A(1,4)=p.a*p.bR*(Sh+(1-p.epsR)*Ph)/p.Nh;
A(2,1)=p.sigma_h; A(2,2)=-rh;
A(3,2)=p.a*p.cR*Sv/p.Nh; A(3,3)=-kv;
A(4,3)=p.sigma_v; A(4,4)=-p.mu_v;

dPhi=A*Phi;
da=[dy;dPhi(:)];
end

function aug=augmented_pulse(aug,p)
y=aug(1:7);
Phi=reshape(aug(8:end),4,4);
y=sensitive_pulse(y,p);
DR=diag([1,1-p.q*p.eR,1,1]);
Phi=DR*Phi;
aug=[y;Phi(:)];
end
