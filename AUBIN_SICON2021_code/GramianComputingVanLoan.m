function [BigMat0,BigMat1] = GramianComputingVanLoan(listT1,listT2,A,B)
% GRAMIANCOMPUTINGVANLOAN outputs two large matrices that are the concatenation of
% the Gramians K(listT1_i,listT2_j) of lisT1*listT2, the output is obtained
% through a "cell2mat" of nn1*nn2 cell_array of N*N matrices.
% BigMat0: (nn1xN)*(nn2xN) for the uncontrolled part K_0 of the kernel
% BigMat1: (nn1xN)*(nn2xN) for the controlled part K_1 of the kernel
% The matrices K(s,t_m) are approximated using Van Loan's trick described in "Computing Integrals Involving the
% Matrix Exponential, CHARLES F. VAN LOAN, TAC,1978". We use the same notations F2 and G2 as in the paper.
% Inputs:
% listT1: 1xnn1, list of time points where the kernel matrices are evaluated
% listT2: 1xnn2, list of time points where the kernel matrices are centered
% A: NxN, B:NxP matrices defining x'=Ax+Bu
% NOTE THAT THE GRAM MATRICES ARE COMPUTED for R=eye(P), Q=0 with the LQR
% classical notations x'Qx+u'Ru.
% The time lists are reordered in the process for Van Loan's formula to hold for the LQR kernel, 
% however the output is put back in the original order.

tic

[listT1,idxSorting1]=sort(listT1);
[listT2,idxSorting2]=sort(listT2);

N=size(A,1); nn1=length(listT1); nn2=length(listT2);
BigMat_cellArr0=cell(nn1,nn1); BigMat_cellArr1=cell(nn1,nn1);
BigMat_cellArrF2_T1=cell(nn1,1); BigMat_cellArrG2_T1=cell(nn1,1);
BigMat_cellArrF2_T2=cell(nn2,1); BigMat_cellArrG2_T2=cell(nn2,1);

if isequal(listT1,listT2) %symmetric case between evaluations and centers, leveraging the Hermitian symmetry of the kernel matrices to reduce the computations.
	for i=1:nn1 %Van Loan's technique for computing Gramians
        temp = expm([A B*B';zeros(N,N) -A']*listT1(i));% REPLACE B*B' BY B*R^{-1}*B' for non diagonal R
        BigMat_cellArrF2_T1{i}=temp(1:N,1:N);
        BigMat_cellArrG2_T1{i}=temp(1:N,N+1:2*N);
    end
    for i=1:nn1 %Computing upperhalf with Van Loan's technique for computing Gramians
		for j=i:nn1
			BigMat_cellArr1{i,j}=BigMat_cellArrG2_T1{i}*BigMat_cellArrF2_T1{j}';
			BigMat_cellArr0{i,j}=BigMat_cellArrF2_T1{i}*BigMat_cellArrF2_T1{j}';
		end
	end
	for i=1:nn1 %Computing lowerhalf by symmetry
		for j=1:i
			BigMat_cellArr1{i,j}=BigMat_cellArr1{j,i}';
			BigMat_cellArr0{i,j}=BigMat_cellArr0{j,i}';
		end
	end
else
    for i=1:nn1 %Van Loan's technique for computing Gramians
        temp = expm([A B*B';zeros(N,N) -A']*listT1(i));
        BigMat_cellArrF2_T1{i}=temp(1:N,1:N);
        BigMat_cellArrG2_T1{i}=temp(1:N,N+1:2*N);
    end
    for i=1:nn2 %Van Loan's technique for computing Gramians
        temp = expm([A B*B';zeros(N,N) -A']*listT2(i));
        BigMat_cellArrF2_T2{i}=temp(1:N,1:N);
        BigMat_cellArrG2_T2{i}=temp(1:N,N+1:2*N);
    end
	for i=1:nn1 %Computing both lower and upper halves with Van Loan's technique for computing Gramians
		for j=1:nn2
			BigMat_cellArr0{i,j}=BigMat_cellArrF2_T1{i}*BigMat_cellArrF2_T2{j}';
			if listT1(i)<= listT2(j)
				BigMat_cellArr1{i,j}=BigMat_cellArrG2_T1{i}*BigMat_cellArrF2_T2{j}';
			else
				BigMat_cellArr1{i,j}=BigMat_cellArrF2_T1{i}*BigMat_cellArrG2_T2{j}';
			end
		end
	end
end
elapsedTime=toc;
disp(['Van Loan: finished Gram ' num2str(elapsedTime) 's']);

idxSorting1_rev(idxSorting1) = 1:length(idxSorting1);
idxSorting2_rev(idxSorting2) = 1:length(idxSorting2);
BigMat0 = cell2mat(BigMat_cellArr0(idxSorting1_rev,idxSorting2_rev));
BigMat1 = cell2mat(BigMat_cellArr1(idxSorting1_rev,idxSorting2_rev));

end

