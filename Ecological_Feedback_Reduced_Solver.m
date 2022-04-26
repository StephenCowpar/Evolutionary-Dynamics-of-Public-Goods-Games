clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 8;
r = 2.76;
d = 0.5;
tmax = 210;
x0 = .8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tspan = [0 tmax];
y0 = [1 x0];
[t,y] = ode45(@(t,y) odefcn(t,y,N,r,d), tspan, y0);

[~,idx1]=min(abs(t-180));
[~,idx2]=min(abs(t-190));
[~,idx3]=min(abs(t-195));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
PS = PLOT_STANDARDS();
figure(1);
fig1_comps.fig = gcf;
hold on
fig1_comps.p1 = plot(t,y(:,1));
fig1_comps.p2 = plot(t,y(:,2));
fig1_comps.p3 = xline(t(idx1));
fig1_comps.p4 = xline(t(idx2));
fig1_comps.p5 = xline(t(idx3));
axis([0 tmax 0 1.2*max(y(:,1))])
hold off
title('');
xlabel('time');
legend([fig1_comps.p1, fig1_comps.p2], '$$H(t)$$', '$$x(t)$$');
legendX = .82; legendY = .87; legendWidth = 0.005; legendHeight = 0.005;
fig1_comps.legendPosition = [legendX, legendY, legendWidth, legendHeight];
set(fig1_comps.p1, 'LineStyle', '-', 'LineWidth', 2.5, 'Color', PS.Blue3);
set(fig1_comps.p2, 'LineStyle', '-', 'LineWidth', 2.5, 'Color', PS.Red3);
set(fig1_comps.p3, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', PS.Orange1);
set(fig1_comps.p4, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', PS.Orange3);
set(fig1_comps.p5, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', PS.Orange5);
STANDARDIZE_FIGURE(fig1_comps);
 
SAVE_MY_FIGURE(fig1_comps, sprintf('testN_%0.0f_r_%0.2f_d_%0.2f_%0.1f_B45_TS.png',N,r,d,x0), 'big');
%% 


PS = PLOT_STANDARDS();
figure(2);
fig2_comps.fig = gcf;
hold on
fig2_comps.p1 = plot(y(:,2),y(:,1));
fig2_comps.p2 = plot(y(1,2),y(1,1),'y.', 'MarkerSize',15);
fig2_comps.p3 = plot(y(end,2),y(end,1),'y.', 'MarkerSize',15);
fig2_comps.p4 = plot(y(idx1,2),y(idx1,1),'y.', 'MarkerSize',15);
fig2_comps.p5 = plot(y(idx2,2),y(idx2,1),'y.', 'MarkerSize',15);
fig2_comps.p6 = plot(y(idx3,2),y(idx3,1),'y.', 'MarkerSize',15);
axis([0 1 0 1.2*max(y(:,1))])
hold off
title('');
xlabel('$$x(t)$$');
ylabel('$$H(t)$$');
%legend([fig1_comps.p1, fig1_comps.p2], '$$H(t)$$', '$$x(t)$$');
%legendX = .82; legendY = .87; legendWidth = 0.005; legendHeight = 0.005;
%fig1_comps.legendPosition = [legendX, legendY, legendWidth, legendHeight];
set(fig2_comps.p1, 'LineStyle', '-', 'LineWidth', 2.5, 'Color', PS.Blue3);
set(fig2_comps.p2,  'MarkerSize', 15, 'Color', PS.Green3);
set(fig2_comps.p3,  'MarkerSize', 15, 'Color', PS.Red3);
set(fig2_comps.p4,  'MarkerSize', 15, 'Color', PS.Orange1);
set(fig2_comps.p5,  'MarkerSize', 15, 'Color', PS.Orange3);
set(fig2_comps.p6,  'MarkerSize', 15, 'Color', PS.Orange5);
STANDARDIZE_FIGURE(fig2_comps);
SAVE_MY_FIGURE(fig2_comps, sprintf('testN_%0.0f_r_%0.2f_d_%0.2f_%0.1f_B45_PS.png',N,r,d,x0), 'big');
%%

PS = PLOT_STANDARDS();
figure(3)
fig3_comps.fig = gcf;
blippy = distsol(y);
s  = 0:.001:1;
hold on
fig3_comps.p1 = plot(s,blippy(1,:));
fig3_comps.p2 = plot(s,blippy(idx1,:));
fig3_comps.p3 = plot(s,blippy(idx2,:));
fig3_comps.p4 = plot(s,blippy(idx3,:));
fig3_comps.p5 = plot(s,blippy(end-1,:));
axis([0 1 0 1])
hold off
xlabel('$$s$$');
ylabel('$$P(s,t)$$');
set(fig3_comps.p1, 'LineStyle', '-', 'LineWidth', 2.5, 'Color', PS.Green3);
set(fig3_comps.p2, 'LineStyle', '-', 'LineWidth', 2.5, 'Color', PS.Orange1);
set(fig3_comps.p3, 'LineStyle', '-', 'LineWidth', 2.5, 'Color', PS.Orange3);
set(fig3_comps.p4, 'LineStyle', '-', 'LineWidth', 2.5, 'Color', PS.Orange5);
set(fig3_comps.p5, 'LineStyle', '-', 'LineWidth', 2.5, 'Color', PS.Red3);
STANDARDIZE_FIGURE(fig3_comps);
SAVE_MY_FIGURE(fig3_comps, sprintf('testN_%0.0f_r_%0.2f_d_%0.2f_%0.1f_B45.png',N,r,d,x0), 'big');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydt = odefcn(~,y,N,r,d)
    dydt = zeros(2,1);
    dydt(1) = y(1) *W(y(2),N,r)*y(2);
    dydt(2) = -(1-y(2))*(Z(y(1))*y(2)*(r-1)*(1-y(2)^(N-1))-d);
end

function val = W(x,N,r)
    val = (-1-(-1+r)*x^(N-1)+(r*(1-x^N))/(N*(1-x))); 
end

function val = Z(x)
    s  = 0:.001:1;
    %G0 = ones(1,length(s));
    %G0 = -(3/5)*(-1-4.*s+4.*x.^2);
     a = 4; b = 5; %mean at a/(a+b).
     G0 = betapdf(s,a,b);
%      pd = makedist('Normal','mu',.1,'sigma',.1);
%      trun = truncate(pd,0,1);
%      G0 = pdf(trun,s);

    val = (trapz(s,s.*G0.*x.^s))/(trapz(s,G0.*x.^s));
end
  
function distret = distsol(y)
    s  = 0:.001:1;
    %  G0 = ones(1,length(s));
    %G0 = -(3/5)*(-1-4.*s+4.*x.^2);
    a = 4; b = 5; %mean at a/(a+b).
         G0 = betapdf(s,a,b);
%     pd = makedist('Normal','mu',.1,'sigma',.1);
%     trun = truncate(pd,0,1);
%     G0 = pdf(trun,s);
    distret = (1-(y(:,2))).*(G0.*(y(:,1)).^s)./trapz(s,G0.*(y(:,1)).^s,2);
end
