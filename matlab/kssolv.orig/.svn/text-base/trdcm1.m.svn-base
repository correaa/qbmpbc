function [Etotvec, X,vtot,rho] = trdcm1(mol,options)
%
%  Usage: [Etotvec, X, vtot, rho] = trdcm1(mol);
%
%         [Etotvec, X, vtot, rho] = trdcm1(mol, options); 
%
%  Purpose:
%     Use trust region enabled direct constrained minization (TRDCM1) 
%     iteration to find the ground state minimum total energy and 
%     the corresponding wave functions.
%
%     Trust region is imposed whenever necessary to ensure the total
%     decreases monotonically.
%
%  Input:
%     mol       a Molecule object
%     options   a structure that can be defined by calling the SETKSOPT
%               funciton.
%  Input:
%     mol       a Molecule object
%     options   a structure that can be defined by calling the SETKSOPT
%               funciton. For example, options = setksopt;
%               defines a default options structure. options can be modified
%               by for example,
%               options = setksopt(options,'maxdcmiter',20,'dcmtol',1e-10);
%
%    The following options fields are used by TRDCM1. Check the setksopt 
%    function (by issuing 'help setksopt' to see what values
%    can be set for these fields.
%      
%    option fields used in TRDCM1:
%
%   verbose      Turn on/off print statements.
%   maxdcmiter   Maximum number of DCM iterations allowed
%   maxcgiter    Maximum number of LOBPCG iterations allowed for solving
%                each linear eigenvalue problem
%   dcmtol       The DCM convergence tolerance. 
%   cgtol        The LOBPCG tolerance.
%   X0           Initial guess of the wavefunctions
%   rho0         Initial guess of the charge density. If not used, atomic
%                charge are used to construct the initial Hamiltonian.
%
%  Output:
%   Etotvec     an array of total energy values computed at each DCM iteration
%   X           final approximate wavefunctions (Wavefun object)
%   vtot        final total (local) potential (3D array)
%   rho         final total charge density (3D array)

%
% check for options
%
if (nargin < 2)
   %
   % use default options
   %
   options = setksopt;
else
   if (~isstruct(options))
      error('The second arguement must be an option structure');
      return;
   end;
end;
%
% parse options 
%
if (any(strcmp(options.verbose,{'off';'OFF';'Off'})))
   verbose = 0;
else
   verbose = 1;
end;
maxdcmiter = options.maxdcmiter;  
maxinerscf = options.maxinerscf;  
maxcgiter  = options.maxcgiter;
dcmtol     = options.dcmtol;
cgtol      = options.cgtol;
X          = options.X0;
rho        = options.rho0;
maxtry     = 10;
fudge      = 1e-12;
fftw('planner','estimate');
%
% fixed trust region parameter
%
sigma = 0.1;
%
% initialize other variables
%
vtot       = [];
Etotvec    = [];
%
% get some size information 
%
n1 = get(mol,'n1');
n2 = get(mol,'n2');
n3 = get(mol,'n3');
n123 = n1*n2*n3;
%
% start timing
%
tstart  = cputime;
%
% figure out the number of occupied state based on whether electron spin is 
% taken into account
%
nspin = get(mol,'nspin');
nocc  = getnel(mol)/2*nspin;
%
% Use the given rho to construction initial Hamiltonian if available
% otherwise always use the atomic charge for the initial construcition
%
if (isempty(rho))
   initrho = 0;
   H = Ham(mol); 
else
   initrho = 1;
   H = Ham(mol,rho); 
end;
%
% generate initial guess of the wavefunction if it is not provided
%
if ( isempty(X) )
   X = genX0(mol,nocc);
   initX = 0;
else
   %
   % check the dimension of X and make sure it is OK
   %
   initX = 1; 
   if ( get(X,'ncols') < nocc )
      fprintf('Error: The number of columns in X is less than nocc\n');
      fprintf('Error: size(X,2) = %d, nocc = %d\n', get(X,'ncols'), nocc);
      return;
   end;
   %
   m1 = get(X,'n1');
   m2 = get(X,'n2');
   m3 = get(X,'n3');
   if (m1 ~= n1 | m2 ~= n2 | m3 ~= n3)
      error('trdcm1: dimension of the molecule does not match that of the wavefunction')
   end;
end;
nxcols = get(X,'ncols');
%
vol = get(mol,'vol');
gmask = FreqMask(mol);
%
iterdcm = 1;
%
%
% extract the ionic and external potential for possible reuse,
% save a copy of total potential for self-consistency check 
%
vion = get(H,'vion');
vin  = get(H,'vtot');
vext = get(H,'vext'); 
%
% calculate Ewald and Ealphat (one time calculation)
%
Ewald     = getewald(mol);
Ealphat   = getealphat(mol);
%fprintf('Ewald   = %11.4e\n', Ewald);
%fprintf('Ealphat = %11.4e\n', Ealphat);
%
% construct a preconditioner for LOBCG
%
prec = genprec(mol);
%
if (~initX || ~initrho)
   %
   % run a few LOBPCG iterations to get a better set of wavefunctions
   %
   [X, ev, lvec, rvec] = lobpcg(H, X, prec, cgtol, maxcgiter,verbose);
   %
   % compute charge density
   %
   rho = getcharge(mol,X);
   %
   % various energy calculation
   %
   Ekin = (2/nspin)*sum(ev(1:nocc));
   Ecor = getEcor(mol,rho,vin,vion,vext);
   %
   % Compute Hartree and exchange correlation energy and potential
   % using the new charge density; update the total potential
   %
   [vhart,vxc,uxc2,rho]=getvhxc(mol,rho);
   vout = getvtot(mol, vion, vext, vhart, vxc);
   %
   % Calculate the potential energy based on the new potential
   %
   Ecoul = getEcoul(mol,abs(rho),vhart);
   Exc   = getExc(mol,abs(rho),uxc2);
   %
   %fprintf('Ekin  = %22.12e\n', Ekin);
   %fprintf('Ecor  = %22.12e\n', Ecor);
   %fprintf('Ecoul = %22.12e\n\n', Ecoul);
   %fprintf('Exc   = %22.12e\n\n', Exc);
   %
   Etot = Ewald + Ealphat + Ekin + Ecor + Ecoul + Exc;
   Etotvec = [Etotvec Etot]; 
   %
   H = set(H,'vtot',vout);  
else
   %
   % Calculate total energy using the given rho and X
   %
   % 1. Calculate the kinetic energy and ionic potential energy 
   %
   KX = applyKIEP(H,X); 
   T = X'*KX; 
   clear KX;
   Ekin = (2/nspin)*trace(T);
   %
   % 2. Calculate the Hartree and exchange-correlation potential energy
   %
   [vhart,vxc,uxc2,rho]=getvhxc(mol,rho);
   Ecoul = getEcoul(mol,abs(rho),vhart);
   Exc   = getExc(mol,abs(rho),uxc2);
   %
   Etot = Ewald + Ealphat + Ekin + Ecoul + Exc;
end
if (verbose==1)
   fprintf('Total energy   = %20.13e\n', Etot);
end; 
%
fprintf('Beging DCM calculation...\n');
while ( iterdcm <= maxdcmiter )
   fprintf('\nDCM iter %2d:\n', iterdcm);
   %
   % calculate the residual 
   %
   HX = H*X;  % there is a more efficient way to compute this.
   T = X'*HX; T = (T+T')/2;  % symmetrize to avoid complex eigenvalues
   R = HX - X*T;
   for j = 1:nocc
      if (verbose == 1)
         fprintf('resnrm = %11.3e\n', norm(R(:,j)));
      end;
      R(:,j)=prec.*R(:,j);
   end;
   %
   % construct the subspace
   %
   Y = [X R];
   if (iterdcm > 1)
      Y = [Y P];
   end;
   ny = get(Y,'ncols');
   %
   % project the kinetic and ionic potential part of the Hamiltonian
   % into the subspace spanned by Y
   %
   KY = applyKIEP(H,Y);
   T = Y'*KY; 
   clear KY;
   %
   B = Y'*Y; B = (B+B')/2;
   %
   % solve the projected problem
   % 
   G = eye(ny);
   for iterscf = 1:maxinerscf
      vnp0 = vhart+vxc;
      % project the nonlinear potential part of the Hamiltonian 
      VY = applyNP(H,Y);
      A = T + Y'*VY;
      A = (A+A')/2;
      BG = B*G(:,1:nocc);
      C = BG*BG';
      [G,D]=eig(A-sigma*C,B,'chol');
      %format short e
      %d = diag(D)
      Ekin = (2/nspin)*trace( G(:,1:nocc)'*T*G(:,1:nocc) );
      X = Y*G(:,1:nocc);
      rho = getcharge(mol,X);
      %
      %  update the the Hartree and exchange correlation potential 
      %  based on the new charge density rho
      % 
      [vhart,vxc,uxc2,rho]=getvhxc(mol,rho);
      %
      % Calculate the potential energy based on the new potential
      %
      Ecoul = getEcoul(mol,abs(rho),vhart);
      Exc   = getExc(mol,abs(rho),uxc2);
      %
      Etot = Ewald + Ealphat + Ekin + Ecoul + Exc;
      if (verbose==1)
         fprintf('   inner_scf = %d, Etot = %20.13e\n', iterscf, Etot);
      end;
      %
      % Update the nonlinear potential only
      %  
      H = set(H,'vnp',vhart+vxc);
   end; 
   Etotvec = [Etotvec Etot];
   %
   % update the total potential
   %
   vout = getvtot(mol, vion, vext, vhart, vxc);
   H = set(H,'vtot',vout);  
   P = Y(:,nocc+1:ny)*G(nocc+1:ny,1:nocc);
   fprintf('norm(vin-vout) = %11.3e\n', ...
            norm(reshape(vhart+vxc-vnp0,n123,1))/norm(reshape(vnp0,n123,1)));
   if (verbose==1)
      fprintf('------\n');
   else
      fprintf('Total energy   = %20.13e\n', Etot);
   end
   %pause;
   %
   iterdcm = iterdcm + 1;
end;
vtot = vout;
timetot = cputime - tstart;
fprintf('Total time used = %11.3e\n', timetot);
HX = H*X(:,1:nocc);
RES = HX - X(:,1:nocc)*(X(:,1:nocc)'*HX);
fprintf('fnorm(HX-XD)    = %11.3e\n', norm(RES,'fro'));
