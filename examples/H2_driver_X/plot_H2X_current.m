% plot_H_current

if (~exist('magnetic'))
  %cell=hdf5read('current.h5','cell');
  %current=hdf5read('current.h5','data');
  %magnetic=hdf5read('magnetic.h5','data');
  
  cell=hdf5read('current90.h5','cell');
  current=hdf5read('current90.h5','data');
  magnetic=hdf5read('magnetic90.h5','data');
  
  %cell=hdf5read('current30.h5','cell');
  %current=hdf5read('current30.h5','data');
  %magnetic=hdf5read('magnetic30.h5','data');
  
  [Nz,dim,Ny,Nx] = size(current);
  Jx = zeros(Nz,Ny,Nx); Jy = Jx; Jz = Jx;
  Jx(:,:,:) = current(:,1,:,:); Jx=permute(Jx,[3,2,1]);
  Jy(:,:,:) = current(:,2,:,:); Jy=permute(Jy,[3,2,1]);
  Jz(:,:,:) = current(:,3,:,:); Jz=permute(Jz,[3,2,1]);
  Bx = zeros(Nz,Ny,Nx); By = Bx; Bz = Bx;
  Bx(:,:,:) = magnetic(:,1,:,:); Bx=permute(Bx,[3,2,1]);
  By(:,:,:) = magnetic(:,2,:,:); By=permute(By,[3,2,1]);
  Bz(:,:,:) = magnetic(:,3,:,:); Bz=permute(Bz,[3,2,1]);
  
  %wf0=hdf5read('wf0.h5','data');
  %wf0real = zeros(Nz,Ny,Nx); wf0imag = wf0real;
  %wf0real(:,:,:) = wf0(:,1,:,:); wf0real=permute(wf0real,[3,2,1]);
  %wf0imag(:,:,:) = wf0(:,2,:,:); wf0imag=permute(wf0imag,[3,2,1]);
  %den0 = wf0real.*wf0real + wf0imag.*wf0imag;
end

x = [0:Nx-1]/Nx * cell(1,1); 
y = [0:Nx-1]/Nx * cell(2,2); 
z = [0:Nx-1]/Nx * cell(3,3); 
X = ones(Nx,Ny,Nz); for i=1:Nx,for j=1:Ny,for k=1:Nz, X(i,j,k) = x(i); end; end; end;
Y = ones(Nx,Ny,Nz); for i=1:Nx,for j=1:Ny,for k=1:Nz, Y(i,j,k) = y(j); end; end; end;
Z = ones(Nx,Ny,Nz); for i=1:Nx,for j=1:Ny,for k=1:Nz, Z(i,j,k) = z(k); end; end; end;

xid = Nx/2 + 1; %mid-plane where atom is

myJx = zeros(Ny,Nz); myJx(:,:) = Jx(xid,:,:);
myJy = zeros(Ny,Nz); myJy(:,:) = Jy(xid,:,:);
myJz = zeros(Ny,Nz); myJz(:,:) = Jz(xid,:,:);

myBx = zeros(Ny,Nz); myBx(:,:) = Bx(xid,:,:);
myBy = zeros(Ny,Nz); myBy(:,:) = By(xid,:,:);
myBz = zeros(Ny,Nz); myBz(:,:) = Bz(xid,:,:);

figure(1)
clf
subplot(2,1,1);

quiver(myJz,myJy);
axis equal
grid
set(gca,'FontSize',14);
title('Qbox / MPBC Hydrogen');
xlabel('x');
ylabel('y');
zlabel('z');
xlim([-15 15]+Nx/2+1);
ylim([-15 15]+Ny/2+1);

subplot(2,1,2);
[sz sy] = meshgrid(Nz/2+1+[-5:5],Ny/2+1);
%[sz sy] = meshgrid(Nz/2+1,Ny/2+1+[-5:5]);
streamline(stream2(myJz,myJy,sz,sy,[0.001,20000]));
axis equal
xlim([-6 6]+Nx/2+1);
ylim([-6 6]+Ny/2+1);
hold on
plot(sz,sy,'ro');
hold off

figure(2);
subplot(2,1,1);
mesh(y,z,-myBx);
disp(sprintf('max Bx on this plane is: %e',max(max(-myBx))));
subplot(2,1,2);
contour(y,z,-myBx);
axis equal


figure(3);
clf
isosurface(abs(Jx.*Jx+Jy.*Jy+Jz.*Jz),5e-7);
axis equal
xlim([0 Nx]);
ylim([0 Ny]);
zlim([0 Nz]);
grid
%(cannot run) quiver3(X,Y,Z,Jx,Jy,Jz);

figure(4);
h2_bond = 1.398397338;
mu0 = 4*pi/(137.036)^2;
au_to_T = 2.35e5;
Bxonline = zeros(Nx,1); Bxonline=-Bx(:,Ny/2+1,Nz/2+1)*mu0;
xgrid = [0:0.01:max(x)]; ygrid=interp1(x,Bxonline,xgrid,'spline');
xH1 = cell(1,1)/2 - h2_bond/2; BxH1 = interp1(x,Bxonline,xH1,'spline');
xH2 = cell(1,1)/2 + h2_bond/2; BxH2 = interp1(x,Bxonline,xH2,'spline');
%plot(x,Bxonline,'.', xgrid,ygrid,'-', xH1,BxH1,'r+', xH2,BxH2,'r+');
plot(x,Bxonline*au_to_T,'.', xgrid,ygrid*au_to_T,'-', xH1,BxH1*au_to_T,'r+', xH2,BxH2*au_to_T,'r+');
xlabel('x (a.u.)');
ylabel('B_x^{ind}  (Tesla)');
Bext = 2*pi/(cell(2,2)*cell(3,3));
title(sprintf('B_x^{ext} = %g  Tesla',Bext*au_to_T)); 
Bind = BxH1; %Bind = mu0*max(max(-myBx));
sigma = Bind/Bext
t1=text(xH2+1,BxH2*au_to_T,sprintf('B_x^{ind} at nucleus = %g T',Bind*au_to_T));
t2=text(xH2+1,BxH2*au_to_T-0.005,sprintf('NMR shift = B_x^{ind}/B_x^{ext} = %.3e',sigma));
