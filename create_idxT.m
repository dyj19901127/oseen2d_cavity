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

idx = idx0;
idx_T = -ones((n+1)*(n+1),1);
for j = 1:(n+1)
  for k = 1:(n+1)
    if mod(j,2)==1 && mod(k,2)==1
      idx_T(((j-1)*(n+1)+k)) = idx;
      idx = idx+1;
    end
  end
end

% Change the reference temperature to 0 
T_eq = mean(T,2) - Tref;

end

