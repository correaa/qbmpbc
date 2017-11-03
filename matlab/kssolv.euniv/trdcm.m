function [Etotvec, X, vtot, rho] = trdcm(mol,options)
%
%  Usage: [Etotvec, X, vtot, rho] = trdcm(mol);
%
%         [Etotvec, X, vtot, rho] = trdcm(mol, options); 
%
%  Purpose:
%     Use trust region enabled direct constrained minization (TRDCM) 
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
%    The following options fields are used by TRDCM. Check the setksopt 
%    function (by issuing 'help setksopt' to see what values
%    can be set for these fields.
%      
%    option fields used in TRDCM:
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
      error('trdcm: dimension of the molecule does not match that of the wavefunction')
   end;
end;
nxcols = get(X,'ncols');
%
vol = get(mol,'vol');
gmask = FreqMask(mol);
%
iterdcm = 1;
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
   [X, ev, lvec, rvec] = lobpcg(H, X, prec, cgtol, maxcgiter, verbose);
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
   % using the new charge density;
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
end;
if (verbose==1)
   fprintf('Total energy   = %20.13e\n', Etot);
end; 
%
fprintf('Beging DCM calculation...\n');
while ( iterdcm <= maxdcmiter )
   %
   % calculate the residual 
   %
   fprintf('\nDCM iter = %d\n', iterdcm);
   % there is a more efficient way to compute this.
   HX = H*X;  
   T = X'*HX; T = (T+T')/2;  % symmetrize to avoid complex eigenvalues
   R = HX - X*T;
   for j = 1:nocc
      if (verbose==1) 
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
   sigma = 0.0;
   iterscf = 1;
   numtry  = 0;
   Etot0   = Etot; 
   while (iterscf <= maxinerscf)
      vnp0 = vhart+vxc;
      % project the nonlinear potential part of the Hamiltonian 
      VY = applyNP(H,Y);
      A = T + Y'*VY;
      A = (A+A')/2;
      BG = B*G(:,1:nocc);
      C = BG*BG';
      C = (C+C')/2;
      if (sigma ~= 0) 
         %
         % trust region imposed, solve a shifted problem
         %
         [G,D]=eig(A-sigma*C,B,'chol');
      else
         [G,D]=eig(A,B,'chol');
      end
      d = diag(D);
      %
      % update wavefunction and charge density
      %   
      X = Y*G(:,1:nxcols);
      rho = getcharge(mol,X);
      %
      %  update the the Hartree and exchange correlation potential 
      %  based on the new charge density rho
      % 
      [vhart,vxc,uxc2,rho]=getvhxc(mol,rho);
      % 
      % kinetic energy
      % 
      Ekin = (2/nspin)*trace( G(:,1:nocc)'*T*G(:,1:nocc) );
      %
      % Calculate the potential energy based on the new potential
      %
      Ecoul = getEcoul(mol,abs(rho),vhart);
      Exc   = getExc(mol,abs(rho),uxc2);
      %
      Etot = Ewald + Ealphat + Ekin + Ecoul + Exc;
      %
      if (verbose==1)
         fprintf('   inner_scf = %d, Etot = %20.13e\n', iterscf, Etot);
      end;
      %
      if (Etot > Etot0 )
         %
         % total energy increased, impose trust region by adjusting sigma
         %
         gaps   = d(2:ny)-d(1:ny-1); 
         gapmax = max(gaps);
         gap0   = d(nocc+1)-d(nocc);
         while (gap0 < 0.9*gapmax & numtry < maxtry)
            if (verbose==1)
               fprintf('increase sigma...');
            end;
            if (sigma == 0)
               sigma = 2*gapmax;
            else
               sigma = 2.0*sigma;
            end; 
            if (verbose==1) fprintf(' sigma = %11.3e\n', sigma); end;
            [G,D] = eig(A-sigma*C,B,'chol');
            d = diag(D); 
            gaps = d(2:ny)-d(1:ny-1); 
            gapmax = max(gaps);
            gap0   = d(nocc+1)-d(nocc);
            numtry = numtry + 1;
         end;
         %
         % --- check Etot again ---
         %
         while (Etot > Etot0 & ...
                abs(Etot-Etot0)>fudge*abs(Etot0) & ...
                numtry < maxtry)
            X = Y*G(:,1:nxcols);
            rho = getcharge(mol,X);
            %
            [vhart,vxc,uxc2,rho]=getvhxc(mol,rho);
            % 
            % kinetic energy
            % 
            Ekin = (2/nspin)*trace( G(:,1:nocc)'*T*G(:,1:nocc) );
            %
            % the potential energy based on the new potential
            %
            Ecoul = getEcoul(mol,abs(rho),vhart);
            Exc   = getExc(mol,abs(rho),uxc2);
            %
            Etot = Ewald + Ealphat + Ekin + Ecoul + Exc;
            if (verbose==1)
               fprintf('   inner_scf = %d, Etot = %20.13e\n', iterscf, Etot);
            end;
            % 
            % increase sigma again if Etot is still larger
            if (Etot > Etot0)
               if (verbose==1)
                  fprintf('increase sigma II...\n');
               end;
               if (sigma~=0) 
                  sigma = 2*sigma;
               else
                  % have to check!!!!
                  sigma = 1.2*gapmax;
               end;  
               if (verbose==1) fprintf(' sigma = %11.3e\n', sigma); end;
               [G,D] = eig(A-sigma*C,B,'chol');
               numtry = numtry + 1;
            end
         end; %while (Etot>Etot0) 
      end; % if (Etot > Etot0)
      %
      % Update the nonlinear potential only
      %  
      H = set(H,'vnp',vhart+vxc);
      Etot0   = Etot; 
      iterscf = iterscf + 1;
   end;
   Etotvec = [Etotvec real(Etot)]; 
   %
   % update the total potential
   %
   vout = getvtot(mol, vion, vext, vhart, vxc);
   H = set(H,'vtot',vout);  
   P = Y(:,nxcols+1:ny)*G(nxcols+1:ny,1:nxcols);
   fprintf('norm(vin-vout) = %11.3e\n', ...
            norm(reshape(vhart+vxc-vnp0,n123,1))/norm(reshape(vnp0,n123,1)));
   if (verbose==1)
      fprintf('------\n');
   else
      fprintf('Total energy   = %20.13e\n', Etot);
   end
   %
   iterdcm = iterdcm + 1;
   %pause;
end;
vtot = vout;
timetot = cputime - tstart;
fprintf('Total time used = %11.3e\n', timetot);
HX = H*X(:,1:nocc);
RES = HX - X(:,1:nocc)*(X(:,1:nocc)'*HX);
fprintf('fnorm(HX-XD)    = %11.3e\n', norm(RES,'fro'));
fprintf('norm(vin-vout)  = %11.3e\n', ...
        norm(reshape(vhart+vxc-vnp0,n123,1))/norm(reshape(vnp0,n123,1)));
%rho = fftshift(real(rho));
%vtot = fftshift(real(vtot));
