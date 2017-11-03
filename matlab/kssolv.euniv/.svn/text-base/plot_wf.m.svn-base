function plot_wf (X, n, mol, isoval)

if (nargin <= 3)
    isoval = -1;
end

psi = get(X,'psi');
idxnz = get(X,'idxnz');
n1 = get(X,'n1');
n2 = get(X,'n2');
n3 = get(X,'n3');
mpbc = get(mol,'mpbc');

if (mpbc==0)
  psi_XYZ = zeros(n1,n2,n3);
  psi_XYZ(idxnz) = psi{n};
  % real space          
  psi_xyz = ifft(ifft(ifft(psi_XYZ,[],1),[],2),[],3);
  figure(1);
  mesh(abs(psi_XYZ(:,:,1)))
  title('wave function in reciprocal space');
  figure(3)
  mesh(abs(psi_xyz(:,:,n3/2)))
  title('wave function in real space');
  if (isoval > 0)
  figure(4);
  clf
  rho_xyz = abs(psi_xyz).^2; 
  max_rho = max(max(max(rho_xyz)));
  isosurface(rho_xyz,max_rho*isoval);
  axis equal
  xlim([0 n1]); ylim([0 n2]); zlim([0 n3]);
  grid
  end
  
else
  % reciprocal space
  psi_XhatZ = zeros(n1*n2,n3);
  psi_XhatZ(idxnz) = psi{n};
  % intermediate space
  psi_xhatZ = ifft(psi_XhatZ,[],1);
  psi_xYZ = reshape(psi_xhatZ,n1,n2,n3);
  % real space          
  psi_xyz = ifft(ifft(psi_xYZ,[],2),[],3);
  figure(1)
  plot(abs(psi_XhatZ(:,1)))
  title('wave function in reciprocal space');
  figure(2)
  plot(abs(psi_xhatZ(:,1)))
  title('wave function in intermediate space');
  figure(3)
  mesh(abs(psi_xyz(:,:,n3/2)))
  title('wave function in real space');
  
  if (isoval > 0)
  figure(4);
  clf
  rho_xyz = abs(psi_xyz).^2; 
  max_rho = max(max(max(rho_xyz)));
  isosurface(rho_xyz,max_rho*isoval);
  axis equal
  xlim([0 n1]); ylim([0 n2]); zlim([0 n3]);
  grid
  end
end