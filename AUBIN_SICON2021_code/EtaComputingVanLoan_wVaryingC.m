function eta_arr = EtaComputingVanLoan_wVaryingC(...
    t_arr_consIneq, delta_arr, nbPointsApprox, Asys, Bsys, C_cellArr_uniq_Ineq, lambda0, lambda1)
%ETACOMPUTINGVANLOAN_wVaryingC outputs eta_arr (#ineq_constraints x 1), 
% which is a vector of eta parameters 1*#unique_ineq_points. 
% This is achieved by approximating the supremum of the definition of eta_m over the whole interval
% taking the supremum over nbPointsApprox grid samples in the interval [t_arr(m)-delta_arr(m),t_arr(m)+delta_arr(m)]
% The matrices K(s,t_m) are approximated using Van Loan's trick described in "Computing Integrals Involving the
% Matrix Exponential, CHARLES F. VAN LOAN, TAC,1978", then the constraint vectors C are applied.
% In practice very few points nbPointsApprox are required. For refined 
% grids t_arr, two points are often enough, the maximum being often
% attained at the border of [t_arr(m)-delta_arr(m),t_arr(m)+delta_arr(m)].

% Inputs:
% t_arr: 1x#time_points, list of time points where to compute eta_arr
% delta_arr: 1x#time_points, list of "radius" of interval around each time point
% Asys: NxN, Bsys:NxP, matrices defining x'=Ax+Bu 
% C_cellArr_uniq_Ineq: cell array of #constraint_points x 1, each cell contains a matrix NxN_cons_m 
% of N_cons_m constraint vectors Nx1 to be applied at point t_m, i.e.
% c_{m,j}^\top x(t_m), j\le N_cons_m.
% lambda0: factor in front of the uncontrolled part of the kernel (related to the choice of norm)
% lambda1: factor in front of the controlled part of the kernel (related to R^{-1} and the choice of norm)
% NOTE THAT THE GRAM MATRICES ARE COMPUTED for R^{-1}=lambda1*eye(P)

% t_arr and C_cellArr should have the same number of elements, in the same
% order, in the sense that t_arr(i) corresponds to the constraints listed
% in C_cellArr{i}. eta_arr has as many of elements as there are constraint
% vectors considered in C_cellArr.

if length(C_cellArr_uniq_Ineq)~=size(t_arr_consIneq,1)
    error('Dimension mismatch between the inequality constraint sets for eta computation')
end

if nargin<8
    lambda1=1;
end
if nargin<7
    lambda0=1;
    lambda1=1;
end

nc=size(blkdiag(C_cellArr_uniq_Ineq{:}),2); N=size(Asys,1); nbPointsApprox=ceil(nbPointsApprox/2)*2;
eta_arr=zeros(nc,1);%length(t_arr),
% eta_arr=zeros(length(t_arr),nc,nn);
% delta_arr=delta.*ones(nn,1);
tic
arr_count=1;
for count=1:length(t_arr_consIneq)
    nn=nbPointsApprox+1; delta=delta_arr(count);
    listT=t_arr_consIneq(count)+linspace(-delta,delta,nn);
    BigMat_cellArrKtutu=cell(nn,1); BigMat_cellArrKttu=cell(nn,1);
    temptt = expm([Asys Bsys*Bsys';zeros(N,N) -Asys']*t_arr_consIneq(count));%K(t,t)
    for i=1:nn %Van Loan's technique for computing Gramians
        temp = expm([Asys Bsys*Bsys';zeros(N,N) -Asys']*listT(i));
        BigMat_cellArrKtutu{i}=lambda0*temp(1:N,1:N)*temp(1:N,1:N)'+lambda1*temp(1:N,N+1:2*N)*temp(1:N,1:N)';%K(t+u,t+u)
        if i<nn/2
        BigMat_cellArrKttu{i}=lambda0*temp(1:N,1:N)*temptt(1:N,1:N)'+lambda1*temp(1:N,N+1:2*N)*temptt(1:N,1:N)';%K(t+u,t) u<0
        else
        BigMat_cellArrKttu{i}=lambda0*temp(1:N,1:N)*temptt(1:N,1:N)'+lambda1*temp(1:N,1:N)*temptt(1:N,N+1:2*N)'; %K(t+u,t) u>0   
        end
    end
    for j=1:size(C_cellArr_uniq_Ineq{count},2) 
        tempC=C_cellArr_uniq_Ineq{count};
        for i=1:nn
         temp_mat=  BigMat_cellArrKtutu{i} + BigMat_cellArrKtutu{ceil(nn/2)} - 2*BigMat_cellArrKttu{i};
        eta_arr(arr_count)= max(eta_arr(arr_count),sqrt(abs(tempC(:,j)'*temp_mat*tempC(:,j))));%computes eta by applying the constraints
        end
        arr_count=arr_count+1;
    end
    
end
elapsedTime=toc;
disp(['Van Loan: finished eta ' num2str(elapsedTime) 's']);
end

% nb_exp=40; lambda=1E-2; c_1=[0 -1 0]'; c_2=[0 0 1]'; C=[c_1, c_2]; delta_arr=linspace(1E-5,0.10,nb_exp);
% eta_arr_1 = EtaComputingVanLoan(repmat(0.5,1,nb_exp),delta_arr, 2, A, B, C, lambda);
% eta_arr_2 = EtaComputingVanLoan(repmat(0.9,1,nb_exp), delta_arr, 2, A, B, C, lambda);
% 
% figure
% hold on
% plot(delta_arr, eta_arr_1(:,1), 'b',delta_arr, eta_arr_1(:,2), 'r')
% plot(delta_arr, eta_arr_2(:,1), '--b',delta_arr, eta_arr_2(:,2), '--r')