%% ========================================================================
% 9. PUBLICATION-QUALITY FIGURE
% ========================================================================

figure('Color','w', ...
       'Position',[80 80 1450 470]);

tl = tiledlayout(1,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

%% ------------------------------------------------------------------------
% PANEL (a): pooled final-5-year resistance
% -------------------------------------------------------------------------

ax1 = nexttile;

bar(ax1,PooledRes,'grouped');

box(ax1,'on');
grid(ax1,'on');

ylim(ax1,[0 105]);
xlim(ax1,[0.5 nf+0.5]);

xticks(ax1,1:nf);
xticklabels(ax1,{'1.00','0.95','0.90'});

xlabel(ax1,'Relative resistant-strain fitness, $f_R$', ...
    'Interpreter','latex');

ylabel(ax1,'Pooled resistance, final 5 years (\%)', ...
    'Interpreter','latex');

title(ax1,'(a) Long-term established resistance', ...
    'FontWeight','normal');

ax1.FontSize = 10;
ax1.LineWidth = 0.8;
ax1.TickDir = 'out';
ax1.GridAlpha = 0.15;

legend(ax1,{ ...
    '$\epsilon_R=0.50$', ...
    '$\epsilon_R=0.70$', ...
    '$\epsilon_R=0.90$'}, ...
    'Interpreter','latex', ...
    'Location','northeast', ...
    'Box','off');

%% ------------------------------------------------------------------------
% PANEL (b): mean malaria prevalence final 5 years
% -------------------------------------------------------------------------

ax2 = nexttile;

bar(ax2,MeanPrev,'grouped');

box(ax2,'on');
grid(ax2,'on');

ymin = min(MeanPrev(:));
ymax = max(MeanPrev(:));

ylim(ax2,[max(0,ymin-0.8), ymax+0.8]);
xlim(ax2,[0.5 nf+0.5]);

xticks(ax2,1:nf);
xticklabels(ax2,{'1.00','0.95','0.90'});

xlabel(ax2,'Relative resistant-strain fitness, $f_R$', ...
    'Interpreter','latex');

ylabel(ax2,'Mean malaria prevalence, final 5 years (\%)', ...
    'Interpreter','latex');

title(ax2,'(b) Long-term malaria burden', ...
    'FontWeight','normal');

ax2.FontSize = 10;
ax2.LineWidth = 0.8;
ax2.TickDir = 'out';
ax2.GridAlpha = 0.15;

%% ------------------------------------------------------------------------
% PANEL (c): contrasting resistance trajectories
% -------------------------------------------------------------------------

ax3 = nexttile;
hold(ax3,'on');

plot(ax3,traj.A.t,traj.A.pR, ...
    'LineWidth',1.8);

plot(ax3,traj.B.t,traj.B.pR, ...
    'LineWidth',1.8);

plot(ax3,traj.C.t,traj.C.pR, ...
    'LineWidth',1.8);

yline(ax3,50,'--', ...
    '50\% threshold', ...
    'Interpreter','latex', ...
    'LabelHorizontalAlignment','left', ...
    'LabelVerticalAlignment','bottom', ...
    'LineWidth',1.0);

xlim(ax3,[0 NYEARS]);
ylim(ax3,[0 105]);

xlabel(ax3,'Time since SMC initiation (years)');

ylabel(ax3,'Resistant-strain frequency (\%)', ...
    'Interpreter','latex');

title(ax3,'(c) Contrasting resistance trajectories', ...
    'FontWeight','normal');

grid(ax3,'on');
box(ax3,'on');

ax3.FontSize = 10;
ax3.LineWidth = 0.8;
ax3.TickDir = 'out';
ax3.GridAlpha = 0.15;

legend(ax3,{ ...
    '$f_R=1.00,\ \epsilon_R=0.50$', ...
    '$f_R=0.95,\ \epsilon_R=0.50$', ...
    '$f_R=0.95,\ \epsilon_R=0.70$'}, ...
    'Interpreter','latex', ...
    'Location','northwest', ...
    'Box','off', ...
    'FontSize',8.5);

hold(ax3,'off');

%% ------------------------------------------------------------------------
% Overall formatting
% -------------------------------------------------------------------------

set(gcf,'Color','w');

% Optional high-resolution export:
% exportgraphics(gcf,'Figure2.png','Resolution',600);
% exportgraphics(gcf,'Figure2.pdf','ContentType','vector');