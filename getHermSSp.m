function [ Vr, Wr ] = getHermSSp( s, A11 ,A12, A21, B1, b, C, c, E11, lf )
%GETHERMSSP Generates the deflating subspace given sigma for IRKA DAE.
%   The following systems are solved 
%
%      |  s_i*E11-A11   A12  | | v_i |     | B1*b_i |
%      |                     | |     |  =  |        |
%      |      A21        0   | |  z  |     |    0   |
%
%      | s_i*E11'-A11'  A21' | | w_i |     | C1'*c_i |
%      |                     | |     |  =  |         |
%      |      A12'       0   | |  z  |     |    0    |
%
%  Then the following deflating matrices are constructed
%
%   Vr = [v_1 v_2 v_3 ... v_r]
%   Wr = [w_1 w_2 w_3 ... w_r]


lf.pmsg(lf.WARN, '=> getHermSSp.m');

lf.pmsg(lf.PED,'    Find the system size.');
[n1,n2] = size(A12);
n = n1+n2;

lf.pmsg(lf.PED,'    Split A into its submatrices');
% A11 = A(1:n1,1:n1);
% A12 = A(1:n1,n1+1:n);
% A21 = A(n1+1:n,1:n1);

lf.pmsg(lf.PED,'    Calculate size of the reduced model.');
r = length(s);

Vr = zeros(n1,r);
Wr = zeros(n1,r);

k = 1;
while k <= r
  lf.pmsg(lf.PED,'Solving for v and w vectors.');
  v = solveDAEVec(s(k),E11, A11, A12, B1, b(:,k), lf );
  w = solveDAEVec(s(k),E11', A11', A12, C', c(:,k), lf );
  
  if (abs(imag(s(k)))/abs(s(k))) < 1e-08 
    lf.pmsg(lf.PED, '      sigma_%d = %f is real.', k, s(k));
    Vr(:,k) = real(v);
    Wr(:,k) = real(w);
    k = k+1;
  else
    lf.pmsg(lf.PED, '      sigma_%d = %f%+fi is complex.', k, real(s(k)), imag(s(k)));
    Vr(:,k:k+1) = [real(v) imag(v)];
    Wr(:,k:k+1) = [real(w) imag(w)];
    k = k+2;
  end
  
end

lf.pmsg(lf.PED,'    Find orthonormal basis for Vr and Wr.');
[Vr,~,~] = svd(Vr,0);
[Wr,~,~] = svd(Wr,0); 

% Wr = Vr;

lf.pmsg(lf.WARN, '<= getHermSSp.m');



end

