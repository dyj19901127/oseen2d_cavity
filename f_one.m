function [ fx,fy ] = f_one( x )
%F_FCN Summary of this function goes here
%   Detailed explanation goes here


fx = zeros(size(x,1));
fy = zeros(size(x,1));

n = size(x,1);
for k=1:n
  if x(k,1)==1
    fy(k) = 1;
  end
end
 


end

