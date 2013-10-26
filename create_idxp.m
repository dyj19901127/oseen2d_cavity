function [ idx_p ] = create_idxp( n, idx0 )
%CREATE_IDXP Create index file for the pressure
%   Detailed explanation goes here

idx = idx0;
idx_p = zeros((n+1)*(n+1),1);
for j = 1:(n+1)
  for k = 1:(n+1)
    if mod(j,2)==1 && mod(k,2)==1
      idx_p(((j-1)*(n+1)+k)) = idx;
      idx = idx+1;
    end
  end
end

end

