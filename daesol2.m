function [ t, x, y ] = daesol2( A, B, C, D, E, tspan, u_fcn, x0)
%DAESOL Summary of this function goes here
%   Detailed explanation goes here

%%% Compute Time Domain Response

    dt = 0.1;
    t = tspan(1):dt:tspan(2);
%    u = zeros(length(t),1);
    u = u_fcn(t);
    
    
    x            = zeros(size(A,1),length(t));
    x(:,1) = x0;
    
    if issparse(A)
      [L, U, P, Q] = lu(E - dt * A);
    else
      [L, U] = lu(E - dt * A);
    end
    
    for i = 2:length(t)
        b      = E*x(:,i-1) + dt*B*u(i)';
        if issparse(A)
          x(:,i) = Q*(U\(L\(P*b)));
        else
          x(:,i) = U\(L\b);
        end
    end
    y = C*x + D*u;


end

