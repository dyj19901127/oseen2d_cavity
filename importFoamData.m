function [u_adv,v_adv,idx_u,idx_p,T,x,e_conn,t] = importFoamData(n,dataLoc)

  if nargin==0
    error('You must supply the size n of the nxn mesh.');
  elseif nargin==1
    dataLoc='data';
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
  T = zeros((n+1)^2,numdirs);

  for k=1:numdirs
    fn = sprintf('%s/%s/U',dataLoc,dataDirs{k});
    fprintf('Importing %s...',fn);
    fid = fopen(fn);
    textscan(fid,'%s %*[^\n]',20);
    tmp = textscan(fid,'%c%f %f %*[^\n]',n*n);
    u(:,k) = interpFoamData(n, tmp{2}, 1,[0,0,0,0]);
    v(:,k) = interpFoamData(n, tmp{3}, 1,[0,0,0,0]);
    fclose(fid);
    fprintf('Complete\n');
    
    % Import T - temperature data
    fn = sprintf('%s/%s/T',dataLoc,dataDirs{k});
    fprintf('Importing %s...',fn);
    if exist(fn,'file')
      fid = fopen(fn);
      textscan(fid,'%s %*[^\n]',20);
      tmp = textscan(fid,'%f',n*n);
      T(:,k) = interpFoamData(n, tmp{1}, 1,[300.5,300,300,299.5]);
      fclose(fid);
      fprintf('Complete\n');
    else
      fprintf('Missing temperature file %s. temperature not imported.',fn);
    end
    
  end
  fprintf('-------------------------------------------------------------------\n');
  fprintf('| Interpolating the centers to the nodes and adding the boundary. |\n');
  fprintf('-------------------------------------------------------------------\n');

  fprintf('-------------------------------------------------------------------\n');
  fprintf('| Creating the mesh...                                            |\n');

  [x, e_conn] = twod_mesh(0,1,0,1,'quadratic',n+1,n+1);

  fprintf('| Completed.                                                      |\n');
  fprintf('-------------------------------------------------------------------\n');
  
  fprintf('-------------------------------------------------------------------\n');
  fprintf('| Creating the u index...                                            |\n');

  [u_adv,v_adv,idx_u] = create_idxu(u,v,1);

  fprintf('| Completed.                                                      |\n');
  fprintf('-------------------------------------------------------------------\n');
 
  fprintf('-------------------------------------------------------------------\n');
  fprintf('| Creating the p index                                          |\n');

  idx_p = create_idxp(n,max(max(idx_u))+1);

  fprintf('| Completed.                                                      |\n');
  fprintf('-------------------------------------------------------------------\n');
  
end



        