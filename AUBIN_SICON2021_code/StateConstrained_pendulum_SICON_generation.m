function [X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq] = StateConstrained_pendulum_SICON_generation(...
    t_samples_consVel,t_samples_consAcc,tgrid_arr)
% STATECONSTRAINED_PENDULUM_SICON_GENERATION takes the list of sampling
% points for the grid and inequalities (which are problem-specific). THE
% PARAMETERS OF THE CONTROL PROBLEM SHOULD BE SET HERE MANUALLY. To
% each unique time point t_m (wherever it appears in the problem) corresponds
% a list of c_{j,m}^\top x(t_m) appearing in the problem. For each c_{j,m} there
% is a real value y_{j,m}, e.g. in the objective y_{j,m}-c_{j,m}x(t_m) or
% in the constraints y_{j,m}+c_{j,m}^\top x(t_m) \ge || = 0. The time points
% should be written with repetition, and set in a matrix "X_..". To each
% corresponds a vector c_{m} Nx1, stored as a row of "C_..", and a real y_m
% stored as a row of "Y_..". The convention for the order in "..._arr" is:
% 1: objective points, 2: inequality points, 3: equality points

% For the output "..._arr_unique" the convention is:
% 1: solution points, 2: objective points, 3: inequality points,
% 4: equality points, 5: grid points

% In the output, the evaluation points are
% unique (to have a better conditionned matrix [K(x_i,x_j)]_{i,j}). The
% inputs "X,Y_arr_unique" are cell arrays 5x1, "C_arr_unique" is 4x1. 
% Each cell of "Y,C_arr_unique" contains a block diagonal matrix of the
% reals/vectors to compute y_{j,m}+c_{j,m}^\top x(t_m) when calling the function
% "..._solver" (WHERE OTHER PARAMETERS SHOULD BE SET)

% C_cellArr_uniq_Ineq is a cell_array containing for each unique t_m of the
% inequality constraints applied to it. It is useful to compute \eta_{m,j}
% for each c_{j,m}.

% Remark: note that X here designates the time arrays, while Y is the state
% output following the classical ML convention Y=c^\top f(X), rather than 
% x=x(t) as in control problems.

t_init=0; t_inter=1/3; t_fin=1; 
vmin=-3; umax=10; umin=-10;
x_init=.5; x_final=0; xDot_init=0; w_init=0; x_inter=0.5;

% n_constraint_types=2; dim_out=3; dim_in=1;

X_arr=cell(3,1);
Y_arr=cell(3,1);
C_arr=cell(3,1);

% GENERATES OBJECTIVE AND CONSTRAINTS

X_samples=t_fin;
Y_samples=[0];

X_obj=X_samples;
Y_obj=Y_samples;
C_obj=[0; -1; 0]';

n_consVel=length(t_samples_consVel);
X_consVel=t_samples_consVel';
C_consVel=repmat([0 1 0],length(t_samples_consVel),1);
Y_consVel=-vmin*ones(n_consVel,1);

n_consAcc=length(t_samples_consAcc);
X_consAccUp=t_samples_consAcc';
C_consAccUp=repmat([0 0 1],length(t_samples_consAcc),1);
Y_consAccUp=-umin*ones(n_consAcc,1);

X_consAccLow=t_samples_consAcc';
C_consAccLow=repmat([0 0 -1],length(t_samples_consAcc),1);
Y_consAccLow=umax*ones(n_consAcc,1);

X_consEq=[t_init; t_init; t_init; t_inter; t_fin];
C_consEq=[[1 0 0];[0 1 0];[0 0 1];[1 0 0];[1 0 0]];
Y_consEq=-[x_init; xDot_init; w_init; x_inter; x_final];

X_arr{1}=X_obj; Y_arr{1}=Y_obj; C_arr{1}=C_obj;
X_arr{2}=[X_consVel;X_consAccUp;X_consAccLow]; 
Y_arr{2}=[Y_consVel;Y_consAccUp;Y_consAccLow]; 
C_arr{2}=[C_consVel;C_consAccUp;C_consAccLow];% Inequality C^\top f(x_m) + Y_m \ge 0
X_arr{3}=X_consEq; Y_arr{3}=Y_consEq; C_arr{3}=C_consEq;% Equality C^\top f(x_m) + Y_m = 0

X_grid=tgrid_arr'; 
Y_grid=[]; 
    
%AUTOMATED PART - FORMATING TO UNIQUE EVALUATION POINTS
[X_arr_uniq,Y_arr_uniq,C_arr_uniq,C_cellArr_uniq_Ineq] =...
    ConvertingConsProblem_toUniquePoints(X_arr,Y_arr,C_arr,X_grid,Y_grid);
end