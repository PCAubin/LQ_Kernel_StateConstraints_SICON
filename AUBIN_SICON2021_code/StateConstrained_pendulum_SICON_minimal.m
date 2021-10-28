% The code below allow to reproduce a minimal numerical result from the
% article:
% P.-C. Aubin-Frankowski "Linearly-constrained Linear Quadratic Regulator from the viewpoint of kernel methods", SICON 2021
% There are some slight differences in the figures due to numerical
% instabilities, the Gram matrices being badly conditioned, and the
% objective being quadratic (which is unadvised for CVX solvers).
% This partial experiment should take about 4s to run on a normal computer.

% CVX http://cvxr.com/cvx/ SHOULD BE INSTALLED, we used a free academic
% license for MOSEK, the solver used in this code 
% It can be retrieved at: https://www.mosek.com/downloads/

% CODE BELOW TO RUN ONCE SO THAT MOSEK+CVX WOULD WORK
% addpath C:\Users\pierr\Mosek\9.2\toolbox\R2015aom
% cd C:\Users\pierr\Desktop\cvx
% cvx_setup
cd 'C:\Users\pierr\OneDrive\Documents\SICON code LQ kernel'
clear all
% close all

% It is possible to use symbolic functions to compute the kernels, but it
% is much slower than using Van Loan's matrix exponential trick
% f_int = @(tau) expm(-A*tau)*B*B'*expm(-A'*tau);    
% K1= @(s,t) expm(s*A)*(integral(f_int, 0, min(s,t), 'ArrayValued',1))*expm(t*A');  
% K0= @(s,t) expm(s*A)*expm(t*A');
%% DEFINITION OF THE PARAMETERS OF THE PROBLEM
tStart = tic;  
g=10;
Asys=[0 1 0; -g 0 1; 0 0 0]; Bsys=[0;0;1]; N=size(Asys,1);
t_init=0; t_inter=1/3; t_fin=1; 
vmin=-3; umax=10; umin=-10;
x_init=.5; x_final=0; xDot_init=0; w_init=0; x_inter=0.5;
% ATTENTION: THESE QUANTITIES REAPPEAR IN THE DEFINITIONS OF THE FUNCTIONS
% "StateConstrained_pendulum_SICON_...", SO MODIFICATIONS SHOULD BE DONE
% THERE TOO
%% DEFINITION OF THE GRID, OF THE KERNEL MATRICES AND SOLUTION
ngrid= 800; tgrid_arr=linspace(0,t_fin,ngrid);
nsamplesVel= 50; t_samples_consVel=linspace(0,t_fin,nsamplesVel);
nsamplesAcc= 200; t_samples_consAcc=linspace(0,t_fin,nsamplesAcc);
bool_SOC=true;
[X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq] = StateConstrained_pendulum_SICON_generation(...
    t_samples_consVel,t_samples_consAcc,tgrid_arr);
%%
[A_sol,C_sol,X_grid,GXgridX] = StateConstrained_pendulum_SICON_solver(...
    X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq,Asys,Bsys,bool_SOC);
%% PLOT OF THE SOLUTION
C_grid_pos = repmat({sparse([1;0;0])},1,length(X_grid));
C_grid_pos = blkdiag(C_grid_pos{:});    
C_grid_vit = repmat({sparse([0;1;0])},1,length(X_grid));
C_grid_vit = blkdiag(C_grid_vit{:});    
C_grid_couple = repmat({sparse([0;0;1])},1,length(X_grid));
C_grid_couple = blkdiag(C_grid_couple{:});    

pos=C_grid_pos'*GXgridX*C_sol*A_sol;
vit=C_grid_vit'*GXgridX*C_sol*A_sol;
couple=C_grid_couple'*GXgridX*C_sol*A_sol;

figure('Renderer', 'painters', 'Position', [100 100 4*217 4*90])
subplot1 =subplot(131); grid
ylim(subplot1,[-0.7 0.7]);

subplot2 =subplot(132); grid
ylim(subplot2,[-4 6]);

subplot3 = subplot(133); grid
ylim(subplot3,[umin-4 umax+4]);

subplot1 =subplot(131);
hold(subplot1, 'on');
plot(X_grid,pos,'r');xlabel('t'); ylabel('x');
plot([0 t_inter t_fin],[x_init x_inter x_final],'o','Color',[0.6350 0.0780 0.1840],'HandleVisibility','off')
subplot2 =subplot(132);
hold(subplot2, 'on');
plot(X_grid,vit,'r');xlabel('t'); ylabel('x'''); 
plot([0],[xDot_init],'o','Color',[0.6350 0.0780 0.1840],'HandleVisibility','off')
text(t_fin,vit(end),num2str(floor(vit(end)*100)/100),'FontSize',8)
subplot3 = subplot(133);
hold(subplot3, 'on');
plot(X_grid,couple,'r','HandleVisibility','off');xlabel('t'); ylabel('w');
text(t_fin,couple(end),num2str(floor(max(couple)*100)/100),'FontSize',8)
subplot(132);
grid
ax0=axis();
h=area([ax0(3) vmin-ax0(3); ax0(3) vmin-ax0(3)],...
    'XData',[0 t_fin],'BaseValue',ax0(3),'ShowBaseLine',0,'LineStyle','none','HandleVisibility','off');
h(2).FaceColor = [0.7 0.7 0.7];
h(2).FaceAlpha = 0.4;

subplot(133);
grid
ax0=axis();
h=area([ax0(3) umin-ax0(3) umax-umin ax0(4)-umax;...
    ax0(3) umin-ax0(3) umax-umin ax0(4)-umax],...
    'XData',[0 t_fin],'BaseValue',ax0(3),'ShowBaseLine',0,'LineStyle','none','HandleVisibility','off');
h(2).FaceColor = [0.7 0.7 0.7];
h(3).FaceColor = [1 1 1];
h(4).FaceColor = [0.7 0.7 0.7];
h(2).FaceAlpha = 0.4;
h(3).FaceAlpha = 0;
h(4).FaceAlpha = 0.4;
subplot(131);
grid

elapsedTime=toc(tStart);
disp(['Total time to perform a simple experiment: ' num2str(elapsedTime) 's']);