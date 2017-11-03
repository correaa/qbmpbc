function [Etotvec, X, vtot, rho] = scf(mol,options)
%
%  Usage: [Etotvec, X, vtot, rho] = scf(mol);
%
%         [Etotvec, X, vtot, rho] = scf(mol,options); 
%
%  Purpose:
%     Use Self Consistent Field (SCF) iteration to find the ground state 
%     minimum total energy and the corresponding wave functions.
%
%  Input:
%     mol       a Molecule object
%     options   a structure that can be defined by calling the SETKSOPT
%               funciton. For example, options = setksopt;
%               defines a default options structure. options can be modified
%               by for example,
%               options = setksopt(options,'maxscfiter',20,'scftol',1e-10);
%
%    The following options fields are used by SCF. Check the setksopt 
%    function (by issuing 'help setksopt' to see what values
%    can be set for these fields.
%      
%    option fields used in SCF:
%
%   verbose      Turn on/off print statements.
%   maxscfiter   Maximum number of SCF iterations allowed
%   maxcgiter    Maximum number of LOBPCG iterations allowed for solving
%                each linear eigenvalue problem
%   scftol       The SCF convergence tolerance. SCF is terminated
%                when norm(vin-vout)/norm(vin) < scftol
%   cgtol        The LOBPCG tolerance.
%   mixtype      The type of potential mixing used in the SCF iteration
%                The possible choices are 'anderson','broyden','broyden1',
%                'pulay','kerker','pulay+kerker'. The default is set to
%                'anderson'.  'broyden' refers to the bad Broyden update.
%                'broyden1' refers to the good Broyden.
%   mixdim       The maximum number potentials to be mixed.
%   betamix      A damping parameter used in Anderson and Broyden mixing.
%                It should be between 0 and 1. The default is 0.8;
%   brank        The rank of the Broyden update
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

%
% check for options
%
if (nargin < 2)
   %
   % use default options
   %
   options = setksopt;
end;
%
% parse options 
%
if (any(strcmp(options.verbose,{'off';'OFF';'Off'})))
   verbose = 0;
else
   verbose = 1;
end;
maxscfiter = options.maxscfiter;  
maxcgiter  = options.maxcgiter;
scftol     = options.scftol;
cgtol      = options.cgtol;
mixtype    = options.mixtype;
mixdim     = options.mixdim;
betamix    = options.betamix;
brank      = options.brank;
X          = options.X0;
rho        = options.rho0;
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
   H = Ham(mol); 
else
   H = Ham(mol,rho); 
end;
%
% generate initial guess of the wavefunction if it is not provided
%
if ( isempty(X) )
   X = genX0(mol,nocc);
else
   %
   % check the dimension of X and make sure it is OK
   %
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
      error('scf: dimension of the molecule does not match that of the wavefunction')
   end;
end;
%
iterscf = 1;
%
% extract the ionic and external potential for possible reuse,
% save a copy of total potential for self-consistency check 
%
vion = get(H,'vion');
vext = get(H,'vext');
vin  = get(H,'vtot');
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
dfmat  = [];
dvmat  = [];
cdfmat = [];

if (~(any(strcmp(mixtype,{'Off';'OFF';'off';'oFF'}))))
   dfmat = zeros(n123,mixdim);
   dvmat = zeros(n123,mixdim);
end;
%
fprintf('Beging SCF calculation...\n');
while ( iterscf <= maxscfiter )
    %
    [X, ev, lvec, rvec] = lobpcg(H, X, prec, cgtol, maxcgiter,verbose);
    %
    % compute charge density
    % 
    rho = getcharge(mol,X);
    %
    % Kinetic energy and some additional energy terms 
    %
    Ekin = (2/nspin)*sum(ev(1:nocc));
    %
    % ionic and external potential energy was included in Ekin
    % along with incorrect Ecoul and Exc. Need to correct them
    % later;
    %
    Ecor = getEcor(mol, rho, vin, vion, vext);
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
    % check progress towards achieving self-consistency in the potential
    % 
    verr = norm(reshape(vout-vin, n123, 1))/norm(reshape(vin,n123,1)); 
    fprintf('\nSCF iter %2d:\n', iterscf);
    if (verbose == 1) 
       fprintf('norm(vout-vin) = %10.3e\n', verr);
    end
    %
    fprintf('Total energy   = %20.13e\n', Etot);
    if (verr < scftol)
       % converged
       fprintf('SCF convergence reached!\n');
       HX = H*X;
       G = X'*HX; G = (G+G')/2;
       RES = HX-X*G;
       if (verbose==1)
          for j = 1:nocc
             fprintf(' resnrm = %11.3e\n', norm(RES(:,j)));
          end
          fprintf('-------------\n\n');
       end
       break;
    end;
    %
    % perform charge mixing
    % 
    vupdated = 0;
    if (iterscf >= 1)
       [vin,dfmat,dvmat,cdfmat] = potmixing(mol, vin,vout,iterscf,mixtype,...
                                  betamix, dfmat, dvmat, cdfmat, mixdim, ...
                                  brank);
       vupdated = 1;
    end
    if (~vupdated)
       vin = vout; 
    end
    %
    % update the total potential and hence the Hamiltonian
    %
    H = set(H,'vtot',vin);
    %
    % check residual error
    %
    HX = H*X;
    G = X'*HX; G = (G+G')/2;
    RES = HX-X*G;
    if (verbose==1)
       for j = 1:nocc
          fprintf(' resnrm = %11.3e\n', norm(RES(:,j)));
       end
       fprintf('-------------\n\n');
    end
    %
    iterscf = iterscf + 1;
    %pause;
end;
vtot = vin;
timetot = cputime - tstart;
fprintf('Total time used = %11.3e\n', timetot);
fprintf('norm(HX-XD)     = %11.3e\n', norm(RES(:,1:nocc),'fro'));
%rho = fftshift(real(rho));
%vtot = fftshift(real(vtot));
