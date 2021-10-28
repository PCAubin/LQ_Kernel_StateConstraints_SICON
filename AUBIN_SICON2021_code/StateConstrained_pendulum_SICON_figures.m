% The code below allow to reproduce the full numerical results of the
% article:
% P.-C. Aubin-Frankowski "Linearly-constrained Linear Quadratic Regulator from the viewpoint of kernel methods", SICON 2021
% There are some slight differences in the figures due to numerical
% instabilities, the Gram matrices being badly conditioned, and the
% objective being quadratic (which is unadvised for CVX solvers).
% The whole experiment should take about 1 min to run on a normal computer.
% It is divided into three part, each corresponding to a given row of
% Figure 1 of the SICON paper.

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
%% 
% Comparison of SOC constraints (guaranteed eta_w) versus discretized constraints (eta_w = 0) for N_P = 200.
figure('Renderer', 'painters', 'Position', [100 100 4*217 4*90])
subplot1 =subplot(131); grid
ylim(subplot1,[-0.7 0.7]);

subplot2 =subplot(132); grid
ylim(subplot2,[-4 6]);

subplot3 = subplot(133); grid
ylim(subplot3,[umin-4 umax+4]);

ngrid= 800; tgrid_arr=linspace(0,t_fin,ngrid);
nsamplesVel= 200; t_samples_consVel=linspace(0,t_fin,nsamplesVel);
nsamplesAcc= 200; t_samples_consAcc=linspace(0,t_fin,nsamplesAcc);
[X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq] = StateConstrained_pendulum_SICON_generation(...
    t_samples_consVel,t_samples_consAcc,tgrid_arr);

list_bool_SOC=[1 0];
for i=1:length(list_bool_SOC)
     bool_SOC=list_bool_SOC(i);
[A_sol,C_sol,X_grid,GXgridX] = StateConstrained_pendulum_SICON_solver(...
    X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq,Asys,Bsys,bool_SOC);

C_grid_pos = repmat({sparse([1;0;0])},1,length(X_grid));
C_grid_pos = blkdiag(C_grid_pos{:});    
C_grid_vit = repmat({sparse([0;1;0])},1,length(X_grid));
C_grid_vit = blkdiag(C_grid_vit{:});    
C_grid_couple = repmat({sparse([0;0;1])},1,length(X_grid));
C_grid_couple = blkdiag(C_grid_couple{:});

pos=C_grid_pos'*GXgridX*C_sol*A_sol;
vit=C_grid_vit'*GXgridX*C_sol*A_sol;
couple=C_grid_couple'*GXgridX*C_sol*A_sol;

color_list=[0 0 0; 0.8500, 0.3250, 0.0980; 1 0 0]; offset_list=.3*[-1 0 1];
offset_nbr=offset_list(i);

subplot1 =subplot(131);
hold(subplot1, 'on');
plot(X_grid,pos,'Color',color_list(i,:));xlabel('t'); ylabel('x');
plot([0 t_inter t_fin],[x_init x_inter x_final],'o','Color',[0.6350 0.0780 0.1840],'HandleVisibility','off')

subplot2 =subplot(132);
hold(subplot2, 'on');
plot(X_grid,vit,'Color',color_list(i,:));xlabel('t'); ylabel('x'''); 
plot([0],[xDot_init],'o','Color',[0.6350 0.0780 0.1840],'HandleVisibility','off')
text(t_fin,vit(end)+offset_nbr,num2str(floor(vit(end)*100)/100),'FontSize',8,'Color',color_list(i,:))%,'Color','red'color_list

subplot3 = subplot(133);
hold(subplot3, 'on');
plot(X_grid,couple,'Color',color_list(i,:),'HandleVisibility','off');xlabel('t'); ylabel('w'); %,Xgap,g*y1+GramD2tm*A+a0+g*a0/2*Xgap.^2,'+'
text(t_fin,couple(end)+3*offset_nbr,num2str(floor(max(couple)*100)/100),'FontSize',8,'Color',color_list(i,:))%,'Color','red'color_list
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
end
lgd=legend('SOC constraints (N_P=200)','discretized constraints');
lgd.Location='southwest';
%%
% Comparison of SOC constraints for varying NP and guaranteed \eta_w.
figure('Renderer', 'painters', 'Position', [100 100 4*217 4*90])
subplot1 =subplot(131); grid
ylim(subplot1,[-0.7 0.7]);

subplot2 =subplot(132); grid
ylim(subplot2,[-4 6]);

subplot3 = subplot(133); grid
ylim(subplot3,[umin-4 umax+4]);

% DEFINITION OF THE GRID AND OF THE KERNEL MATRICES

list_nb_samplingPts=[200, 400, 800];
for i=1:length(list_nb_samplingPts)
ngrid= 800; tgrid_arr=linspace(0,t_fin,ngrid);
nsamplesVel= list_nb_samplingPts(i); t_samples_consVel=linspace(0,t_fin,nsamplesVel);
nsamplesAcc= list_nb_samplingPts(i); t_samples_consAcc=linspace(0,t_fin,nsamplesAcc);
bool_SOC=true;
[X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq] = StateConstrained_pendulum_SICON_generation(...
    t_samples_consVel,t_samples_consAcc,tgrid_arr);
[A_sol,C_sol,X_grid,GXgridX] = StateConstrained_pendulum_SICON_solver(...
    X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq,Asys,Bsys,bool_SOC);

C_grid_pos = repmat({sparse([1;0;0])},1,length(X_grid));
C_grid_pos = blkdiag(C_grid_pos{:});    
C_grid_vit = repmat({sparse([0;1;0])},1,length(X_grid));
C_grid_vit = blkdiag(C_grid_vit{:});    
C_grid_couple = repmat({sparse([0;0;1])},1,length(X_grid));
C_grid_couple = blkdiag(C_grid_couple{:});

pos=C_grid_pos'*GXgridX*C_sol*A_sol;
vit=C_grid_vit'*GXgridX*C_sol*A_sol;
couple=C_grid_couple'*GXgridX*C_sol*A_sol;

color_list=[0 0 0; 0.8500, 0.3250, 0.0980; 1 0 0]; offset_list=.2*[-1 0 1];
offset_nbr=offset_list(i);

subplot1 =subplot(131);
hold(subplot1, 'on');
plot(X_grid,pos,'Color',color_list(i,:));xlabel('t'); ylabel('x');
plot([0 t_inter t_fin],[x_init x_inter x_final],'o','Color',[0.6350 0.0780 0.1840],'HandleVisibility','off')

subplot2 =subplot(132);
hold(subplot2, 'on');
plot(X_grid,vit,'Color',color_list(i,:));xlabel('t'); ylabel('x'''); 
plot([0],[xDot_init],'o','Color',[0.6350 0.0780 0.1840],'HandleVisibility','off')
text(t_fin,vit(end)+offset_nbr,num2str(floor(vit(end)*100)/100),'FontSize',8,'Color',color_list(i,:))%,'Color','red'color_list

subplot3 = subplot(133);
hold(subplot3, 'on');
plot(X_grid,couple,'Color',color_list(i,:),'HandleVisibility','off');xlabel('t'); ylabel('w'); %,Xgap,g*y1+GramD2tm*A+a0+g*a0/2*Xgap.^2,'+'
text(t_fin,couple(end)+3*offset_nbr,num2str(floor(max(couple)*100)/100),'FontSize',8,'Color',color_list(i,:))%,'Color','red'color_list
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
end
lgd=legend(['N_P=' num2str(list_nb_samplingPts(1))],['N_P=' num2str(list_nb_samplingPts(2))],...
    ['N_P=' num2str(list_nb_samplingPts(3))]);%,'FRecon','LRecon'
lgd.Location='southwest';
%% 
% Comparison of SOC constraints for varying \eta_w and NP = 200.
figure('Renderer', 'painters', 'Position', [100 100 4*217 4*90])
subplot1 =subplot(131); grid
ylim(subplot1,[-0.7 0.7]);

subplot2 =subplot(132); grid
ylim(subplot2,[-4 6]);

subplot3 = subplot(133); grid
ylim(subplot3,[umin-4 umax+4]);

% DEFINITION OF THE GRID AND OF THE KERNEL MATRICES
ngrid= 800; tgrid_arr=linspace(0,t_fin,ngrid);
nsamplesVel= 200; t_samples_consVel=linspace(0,t_fin,nsamplesVel);
nsamplesAcc= 200; t_samples_consAcc=linspace(0,t_fin,nsamplesAcc);
bool_SOC=false;
[X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq] = StateConstrained_pendulum_SICON_generation(...
    t_samples_consVel,t_samples_consAcc,tgrid_arr);
lambda=1E-2;
list_eta=[0.05,1E-2,1E-3]/sqrt(lambda);
for i=1:length(list_eta)
     bool_SOC=true;
C_prescribed=[0 0; 0 0; 1 -1]; eta_prescribed=[list_eta(i),list_eta(i)];
[A_sol,C_sol,X_grid,GXgridX] = StateConstrained_pendulum_SICON_solver(...
    X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq,Asys,Bsys,bool_SOC,C_prescribed,eta_prescribed);

C_grid_pos = repmat({sparse([1;0;0])},1,length(X_grid));
C_grid_pos = blkdiag(C_grid_pos{:});    
C_grid_vit = repmat({sparse([0;1;0])},1,length(X_grid));
C_grid_vit = blkdiag(C_grid_vit{:});    
C_grid_couple = repmat({sparse([0;0;1])},1,length(X_grid));
C_grid_couple = blkdiag(C_grid_couple{:});

pos=C_grid_pos'*GXgridX*C_sol*A_sol;
vit=C_grid_vit'*GXgridX*C_sol*A_sol;
couple=C_grid_couple'*GXgridX*C_sol*A_sol;

color_list=[0 0 0; 0.8500, 0.3250, 0.0980; 1 0 0]; offset_list=.3*[-1 0 1];
offset_nbr=offset_list(i);

subplot1 =subplot(131);
hold(subplot1, 'on');
plot(X_grid,pos,'Color',color_list(i,:));xlabel('t'); ylabel('x');
plot([0 t_inter t_fin],[x_init x_inter x_final],'o','Color',[0.6350 0.0780 0.1840],'HandleVisibility','off')

subplot2 =subplot(132);
hold(subplot2, 'on');
plot(X_grid,vit,'Color',color_list(i,:));xlabel('t'); ylabel('x'''); 
plot([0],[xDot_init],'o','Color',[0.6350 0.0780 0.1840],'HandleVisibility','off')
text(t_fin,vit(end)+offset_nbr,num2str(floor(vit(end)*100)/100),'FontSize',8,'Color',color_list(i,:))%,'Color','red'color_list

subplot3 = subplot(133);
hold(subplot3, 'on');
plot(X_grid,couple,'Color',color_list(i,:),'HandleVisibility','off');xlabel('t'); ylabel('w'); %,Xgap,g*y1+GramD2tm*A+a0+g*a0/2*Xgap.^2,'+'
text(t_fin,couple(end)+3*offset_nbr,num2str(floor(max(couple)*100)/100),'FontSize',8,'Color',color_list(i,:))%,'Color','red'color_list
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
end
lgd=legend(['\eta_w= ' num2str(floor(list_eta(1)*1000)/1000) ' (guaranteed)'],['\eta_w= ' num2str(list_eta(2))],...
    ['\eta_w= ' num2str(list_eta(3))]);
lgd.Location='southwest';

elapsedTime=toc(tStart);
disp(['Total time to perform all SICON experiments: ' num2str(elapsedTime) 's']);


% fig = gcf; fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'Pendulum_SOC_fixedN_varyingEta','-dpdf')