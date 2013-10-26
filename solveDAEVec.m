function [ v ] = solveDAEVec( s, E11, A11, A12, B1, b, lf )
%SOLVEDAEVEC This a a helper function used to determine the 
%   This function solves the following system
%
%      -                 - -   -     -      -
%      | s*E11-A11  A12  | | v |     | B1*b |
%      |                 | |   |  =  |      |
%      |    A12'     0   | | z |     |   0  |
%      -                 - -   -     -      -
%

lf.pmsg(lf.PED, '=> solveDAEVec.m');

[n1,n2] = size(A12);

lf.pmsg(lf.PED, '  System size = n1 + n2 = %d + %d = %d.',n1,n2,n1+n2);


lf.pmsg(lf.PED, '  Creating matrices to solve for v.');
Av = [(s*E11)-A11, A12; A12', zeros(n2)];
bv = [B1*b; zeros(n2,1)];


lf.pmsg(lf.PED, '  Solving for v.');
v = Av\bv;

v = v(1:n1);

lf.pmsg(lf.PED, '<= solveDAEVec.m');

end

