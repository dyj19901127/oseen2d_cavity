function [  ] = showFoamData( x, data, num)
%SHOWFOAMDATA Summary of this function goes here
%   Detailed explanation goes here

n = round(sqrt(length(x)));
X = reshape(x(:,1),n,n);
Y = reshape(x(:,2),n,n);

if num == 0
  num = size(data,2);
end

for k = 1:num
  Z = reshape(data(:,k),n,n);
  figure;
  surf(X,Y,Z);
end

end

