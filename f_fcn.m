function [ fx,fy ] = f_fcn( x )
%F_FCN Summary of this function goes here
%   Detailed explanation goes here


fx = (x(:,2)-0.5);
fy = (0.5-x(:,1));
 


end

