function [T_eq,idx_T] = create_idxT(T, n, idx0, Tref)
%CREATE_IDXU creates U and idx_u for oseen_2d_matrices
%
% Inputs:
%   T = Temperature data
%   n = size of matrix
%   idx0 = last index used.
%
% Ouputs:
%   T_eq = averaged temperature data over the time steps.
%   idx_T = indices for temperature vector
%

idx = idx0:(idx0+n^2-2);

lbnd = -ones(1,n+1);
rbnd = -2*ones(1,n+1);

idx_T = [lbnd;reshape(idx,n-1,n+1);rbnd];

idx_T = idx_T(:);



% Change the reference temperature to 0 
T_eq = mean(T,2) - Tref;

end

