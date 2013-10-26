function [ s ] = set_init_shifts( distr, n, distr_type, lf )
%SET_INIT_SHIFTS This sets the initial shift points for IRKA
%
% This produces a 1xn vector of shift points strictly in the right half
% plane.  If im1,im2 = 0 then the shifts returned are real.  
%   Inputs:
%
%   Outputs:
%
%--------------------------------------------------------------------------
% Author: Alan Lattimer
% email: alattime at vt dot edu
% Date: April 2013
%
if nargin < 4
  error('Not enough input arguments');  
end

re1 = distr(1); re2 = distr(2);
im1 = distr(3); im2 = distr(4);

lf.pmsg(lf.WARN, '=> set_init_shifts.m');
s = ones(1,n);

if lf.lvlck(lf.WARN)
  figure('name','Initial Shifts','NumberTitle','off');
end

if distr_type==1 || (im1==0 && im2==0)
  lf.pmsg(lf.WARN, 'Returning Real sigma.');
  s = logspace(re1, re2, n);
  if lf.lvlck(lf.WARN)
    semilogx(s,zeros(size(s)),'*');
    xlabel('Real');
    ylabel('Imaginary');
  end;
else
  lf.pmsg(lf.ERR, 'Setting the shifts to 1, complex shifts not ready.');
  if mod(n,2)==1
    s(n) = re2;
    n = n-1;
  end
  switch distr_type
    case 2
      z = complex(logspace(re1,re2,n/2),logspace(im1,im2,n/2));
      z = [z;conj(z)];
      s(1:n) = z(:);
    case 3
      x = wrev(logspace(-1,pi,n/2)-0.05);
      y = 0.5*(exp(1i*x)+exp(2*1i*x));plot(y,'.');
      [~,ord_re] = sort(real(y));
      [~,ord_im] = sort(imag(y));
      new_re = logspace(re1,re2,n/2);
      new_im = logspace(im1,im2,n/2);
      z = complex(new_re(ord_re),new_im(ord_im));
      z = [z;conj(z)];
      s(1:n) = z(:);
    case 4
      scale = wrev(logspace(re1,re2,n/2));
      x = logspace(-1,pi,n/2)-0.05;
      y = exp(1i*x);plot(y,'.');
      z = complex((real(y)./2 + 0.5).*scale,imag(y));
      z = [z;conj(z)];
      s(1:n) = z(:);
    case 5
      z = complex(logspace(-2.5,-0.5,n/2),logspace(-0.5,1.5,n/2));
      z = [z;conj(z)];
      s(1:n) = z(:);
    otherwise
      z = complex(logspace(re1,re2,n/2),logspace(im1,im2,n/2));
      z = [z;conj(z)];
      s(1:n) = z(:);
  end
  if lf.lvlck(lf.WARN)
    loglog(s,'*');
    xlabel('Real');
    ylabel('Imaginary');
  end
  
end


s = s';

lf.pmsg(lf.WARN, '<= set_init_shifts.m');

end

