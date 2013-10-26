function [v] = interpFoamData(n,data,add_bnd,bnd)
%AVGFOAMDATA Takes center data from OpenFOAM and interpolates it to the vertices
%   input:
%     n         number of elements in each direction on the unit square  
%     data      OF data as column vector
%     add_bnd   add boundary nodes 
%                   0 - No boundary nodes
%                   1 - Add boundary nodes
%     bnd       [l, t, b, r] vector containing the boundary values for 
%
%                           t
%                     -------------
%                     |           |
%                     |           | 
%                   l |           | r
%                     |           | 
%                     -------------
%                           b
%  
% 



k = 1;

vtmp = zeros(n-1^2,1);


for idx=1:(n-1)^2
  vtmp(idx) = (data(k)+data(k+1)+data(k+n)+data(k+n+1))/4;
  if mod(idx, n-1)==0
    k = k+2;
  else
    k = k+1;
  end
  
end

if add_bnd == 0
  v = vtmp;
else
  v1 = bnd(1)*ones(n+1,1);
  v4 = bnd(4)*ones(n+1,1);
  if add_bnd==1
    v2 = bnd(2)*ones(1,n-1);
    v3 = bnd(3)*ones(1,n-1);
  else
    v3 = ((data(1:n-1)+data(2:n))./2)';
    v2 = ((data(n^2-n+1:n^2-1) + data(n^2-n+2:n^2))./2)';
  end
  
  vtmp = [v1, [v2;flipud(reshape(vtmp,n-1,n-1)');v3], v4];
  v = reshape( flipud(vtmp)',[(n+1)^2,1]);
end

end




