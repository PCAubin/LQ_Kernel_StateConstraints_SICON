function [A_sol,C_sol,X_grid,GXgridX] = StateConstrained_pendulum_SICON_solver(...
    X_arr,Y_arr,C_arr,C_cellArr_uniq_Ineq,Asys,Bsys,bool_SOC,C_prescribed,eta_prescribed)
% STATECONSTRAINED_PENDULUM_SICON_SOLVER computes the coefficients A_sol
% (N_coeff x 1) solving the constrained problem. THE PARAMETERS OF THE
% CONTROL PROBLEM SHOULD BE SET HERE MANUALLY. The related vectors (based
% on Gauss pivot) C_sol (#active_points x N_coeff) and the grid X_grid
% (#grid_pts x N) are passed from the inputs. GXgridX (#grid_pts x
% active_points) is the Gram matrix between the test points and the active
% points. The main difference with
% STATECONSTRAINED_PENDULUM_SICON_GENERATION is that here Asys,Bsys are
% taken as inputs to compute the LQ kernel Gram matrices. Hence one should
% also set the ratios 

% Inputs:
% X_arr,Y_arr,C_arr,C_cellArr_uniq_Ineq: cell array of block diagonal 
% matrices obtained as output of STATECONSTRAINED_PENDULUM_SICON_GENERATION
% Asys: NxN, Bsys:NxP, matrices defining x'=Ax+Bu
% bool_SOC: activates the SOC constraints in the problem to solve
% C_prescribed: Nx#prescribed_constraints, matrix of constraint for which
% the value of eta should be found in the vector eta_prescribed
% (#prescribed_constraints x 1) instead of its official value as returned
% by ETACOMPUTINGVANLOAN_wVaryingC.

% YOU SHOULD SET MANUALLY THE RATIOS BETWEEN K0, K1, and the objective, e.g.
% if one minimizes -xDot(T)*lambda_obj + lambdaU \int|u|^2 + lambda_init |x(0)|^2
% then lambda0=1/sqrt(lambda_init), lambda1=1/sqrt(lambdaU) (corresponding to
% a diagonal weight R^{-1} in the kernel matrices).

lambda0=1; lambda1=1E2; lambda_obj=1E2;

% Remark: note that X here designates the time arrays, while Y is the state
% output following the classical ML convention Y=c^\top f(X), rather than 
% x=x(t) as in control problems.

X_sol=X_arr{1}; C_sol=C_arr{1};
X_samples=X_arr{2}; Y_samples=Y_arr{2}; C_samples=C_arr{2};% objective
X_consIneq=X_arr{3}; Y_consIneq=Y_arr{3}; C_consIneq=C_arr{3};% inequality constraints
X_consEq=X_arr{4}; Y_consEq=Y_arr{4}; C_consEq=C_arr{4};% equality constraints
X_grid=X_arr{end}; Y_grid=Y_arr{end};% grid definition

if nargin<8
    C_prescribed=[];
    eta_prescribed=[];
end
if nargin<9
    C_prescribed=[];
    eta_prescribed=[];
end

[GXX0,GXX1]=GramianComputingVanLoan(reshape(X_sol,1,[]),...
    reshape(X_sol,1,[]),Asys,Bsys); %computed for R=Id

GXX=lambda0*GXX0+lambda1*GXX1;
tol=1E-8; sqrtGXX=chol(GXX+tol*eye(size(GXX)));

[X_consIneq_sort,idx_X_consIneq_sort]=sort(X_consIneq);
delta_arr=max(abs([diff(X_consIneq_sort);0]),abs([0;diff(flip(X_consIneq_sort))]))/2;
idx_X_consIneq_sort_rev(idx_X_consIneq_sort) = 1:length(idx_X_consIneq_sort);

if bool_SOC
    nbPointsApprox=2;
eta_arr = EtaComputingVanLoan_wVaryingC(X_consIneq, delta_arr(idx_X_consIneq_sort_rev),...
    nbPointsApprox,Asys, Bsys, C_cellArr_uniq_Ineq, lambda0, lambda1);
    for j=1:size(C_prescribed,2)
        eta_arr(ismember(cell2mat(C_cellArr_uniq_Ineq')',C_prescribed(:,j)','rows'))...
            =eta_prescribed(j);
    end
end
mat_obj=C_samples'*GXX*C_sol;
mat_consIneq=C_consIneq'*GXX*C_sol;
mat_consEq=C_consEq'*GXX*C_sol;
mat_consNorm=sqrtGXX*C_sol;

mat_objNorm=C_sol'*GXX1*C_sol;

n_coeff=size(C_sol,2); tol=1E-8;
cvx_begin
cvx_solver mosek
    variables A(n_coeff,1)
    minimize((full(sum(Y_samples))'+mat_obj*A)*lambda_obj...
    +quad_form(A, mat_objNorm+tol*eye(size(mat_objNorm))))%
    subject to %
        0 == full(sum(Y_consEq))'+ mat_consEq*A;
       if bool_SOC
       eta_arr*norm(mat_consNorm*A)<= full(sum(Y_consIneq))'+ mat_consIneq*A;
       else
       0 <= full(sum(Y_consIneq))'+ mat_consIneq*A;
       end
cvx_end
A_sol=A;
[GXgridX0,GXgridX1]=GramianComputingVanLoan(reshape(X_grid,1,[]),...
    reshape(X_sol,1,[]),Asys,Bsys);
GXgridX=lambda0*GXgridX0+lambda1*GXgridX1;
end

