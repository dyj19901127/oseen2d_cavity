function [u_adv,v_adv,T_eq,idx_u,idx_T,idx_p,x,e_conn,t] = importFoamData(n,dataLoc,Tref,lf)

  if nargin==0
    error('You must supply the size n of the nxn mesh.');
  elseif nargin==1
    dataLoc='data';
  elseif nargin==2
    skipT = 1;
    Tref = 0;
  else
    skipT = 0;
  end

%   h = 1/n;
% 
%   x = zeros(n^2,2);
% 
%   for j=1:n
%     for k=1:n
%       x(n*(j-1) + k,1) = h*(k-0.5);
%       x(n*(j-1) + k,2) = h*(j-0.5);
%     end
%   end

  dataDirs = ls(dataLoc);
  sort(dataDirs);

  dataDirs = textscan(dataDirs,'%s');
  dataDirs = sort(dataDirs{1});
  t = cellfun(@str2num,dataDirs);

  numdirs = size(dataDirs,1);

  u = zeros((n+1)^2,numdirs);
  v = zeros((n+1)^2,numdirs);
  if skipT
    T = 0;
  else
    T = zeros((n+1)^2,numdirs);
  end

  lf.pmsg(lf.ERR,'Loading OpenFOAM data from...')
  lf.pmsg(lf.ERR,' %s',dataLoc)
  for k=1:numdirs
    fn = sprintf('%s/%s/U',dataLoc,dataDirs{k});
    lf.pmsg(lf.WARN,'++ Importing %s...',fn);
    fid = fopen(fn);
    textscan(fid,'%s %*[^\n]',20);
    tmp = textscan(fid,'%c%f %f %*[^\n]',n*n);
    u(:,k) = interpFoamData(n, tmp{2}, 1,[0,0,0,0]);
    v(:,k) = interpFoamData(n, tmp{3}, 1,[0,0,0,0]);
    fclose(fid);
    lf.pmsg(lf.WARN,'   Complete.');
    
    if ~skipT
      % Import T - temperature data
      fn = sprintf('%s/%s/T',dataLoc,dataDirs{k});
      lf.pmsg(lf.WARN,'++ Importing %s...',fn);
      if exist(fn,'file')
        fid = fopen(fn);
        textscan(fid,'%s %*[^\n]',20);
        tmp = textscan(fid,'%f',n*n);
        T(:,k) = interpFoamData(n, tmp{1}, 2,[Tref+0.5, Tref, Tref, Tref-0.5]);
        fclose(fid);
        lf.pmsg(lf.WARN,'   Complete.');
      else
        lf.pmsg(lf.ERR,'Missing temperature file %s. temperature not imported.',fn);
      end
    end
    
  end
  lf.pmsg(lf.ERR,'Interpolating the centers to the nodes and adding the boundary.');

  lf.pmsg(lf.ERR,'+ Creating the mesh...');
  [x, e_conn] = twod_mesh(0,1,0,1,'quadratic',n+1,n+1);
  lf.pmsg(lf.ERR,'+ Completed.');

  lf.pmsg(lf.ERR,'+ Creating the u index...');
  [u_adv,v_adv,idx_u] = create_idxu(u,v,1);
  lf.pmsg(lf.ERR,'+ Completed.');

  if skipT==1
    lf.pmsg(lf.ERR,'+ Skipping T idx build.');
    T_eq = T;
    idx_T = max(max(idx_u));
  else
    lf.pmsg(lf.ERR,'+ Creating the T_eq vector and T index');
    [T_eq,idx_T] = create_idxT(T,n,max(max(idx_u))+1,Tref);
    lf.pmsg(lf.ERR,'+ Completed.');
  end
  


  lf.pmsg(lf.ERR,'+ Creating the p index');
  idx_p = create_idxp(n,max(max(idx_T))+1);

  lf.pmsg(lf.ERR,'+ Completed.');
  
end



        