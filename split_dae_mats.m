function [E11, A11, A12, A21, B1, B2, C1, C2] = split_dae_mats(A, B, C, E, nv)
%SPLIT_DAE_MATS Splits Matrices into DAE submatrices
%   Usage:
%     [E11, A11, A12, A21, B1, B2, C1, C2] = split_dae_mats(A, B, C, E, nv)
%
%

E11 = E(1:nv,1:nv);
A11 = A(1:nv,1:nv);
A12 = A(1:nv,(nv+1):end);
A21 = A((nv+1):end,1:nv);
B1 = B(1:nv,:);
B2 = B((nv+1):end,:);
C1 = C(:,1:nv);
C2 = C(:,(nv+1):end);


end

