function [ u ] = u_fcn( t )
%U_FCN Control function
%   This is a control function for the free convection cavity.

u = ones(size(t));
u(t>3&t<=4) = -10;
u(t>4&t<=6) = 5;
u(t>6) = -0.9;


end

