function H = Ham(varargin)
%
% Hamiltonian constructor
%
% usage:  H = Ham();         constructs a null Hamiltonian
%         H = Ham(mol);      constructs an initial Hamiltonian for Molecule mol
%                            This will use the sum of atomic charges to 
%                            approximate the total density
%
%         H = Ham(mol,rho);  constructs an initial Hamiltonian from
%                            the Molecule mol and precalculated charge
%                            density rho.
%
switch (nargin)
  case 0
     H.gkk  = [];
     H.gm   = [];
     H.vtot = [];
     H.vion = [];
     H.vext = []; 
     H.vnp  = [];
     H.wqmat = []; 
     H.wqsign = [];
     H.mpbc = [];
     H.mpbc_denom = [];
     H.supercell = [];
     H = class(H, 'Ham');
  case 1
     %
     % the input arguement is a Molecule object
     %
     H.gkk  = [];
     H.gm   = [];
     H.vtot = [];
     H.vion = [];
     H.vext = []; 
     H.vnp  = [];
     H.wqmat = []; 
     H.wqsign = []; 
     H.mpbc = [];
     H.mpbc_denom = [];
     H.supercell = [];
     mol = varargin{1}; 
     if ( isa(mol,'Molecule') )
        gmask  = FreqMask(mol);
        H.gkk  = get(gmask,'gkk');
        H.gm   = get(gmask,'gm');
        pseudovar = pseudoinit(mol);
        atomlist = get(mol,'atomlist');
        %
        % compute the local ionic potential
        % 
        [H.vion,vionT,rho] = getvion(mol,pseudovar);
        %
        % compute nonlocal ionic potential
        %
        [H.wqmat,H.wqsign] = getwq(mol,pseudovar);
        clear pseudovar;
        % 
        % Calculate Hartree and exchange-correlation potential
        % from the known density rho
        %
        [vhart,vxc,uxc2,rho]=getvhxc(mol,abs(rho));
        %
        % Save a copy of the Hartree and exchange-correlation
        % separately for DCM updating purposes
        %
        H.vnp = vhart+vxc;
        %
        % there may be other external potential
        % 
        H.vext = get(mol,'vext'); 
        %
        % Sum up all the local potential (sinc interpolation)
        % Note that vtot does not include non-local ionici potential
        % 
        H.vtot = getvtot(mol, H.vion, H.vext, vhart, vxc);
        H.rho  = rho;
        H.mpbc = get(mol,'mpbc');
        H.mpbc_denom = get(mol,'mpbc_denom');
        H.supercell = get(mol,'supercell');
        H = class(H,'Ham'); 
     else 
        error('The input must be a Molecule object');
     end; 
  case 2
     %
     % the first input arguement is mol, the second is rho
     %
     H.gkk  = [];
     H.gm   = [];
     H.vtot = [];
     H.vion = [];
     H.vext = []; % bhlee
     H.vnp  = [];
     H.wqmat = []; 
     H.wqsign = []; 
     H.mpbc = [];
     H.mpbc_denom = [];
     H.supercell = [];
     mol = varargin{1};  
     if ( isa(mol,'Molecule') )
        n1     = get(mol,'n1');
        n2     = get(mol,'n2');
        n3     = get(mol,'n3');
        gmask  = FreqMask(mol);
        H.gkk  = get(gmask,'gkk');
        H.gm   = get(gmask,'gm');
        rho = varargin{2};
        if ( isnumeric(rho) )
           [m1,m2,m3] = size(rho);
           if ( m1 == n1 & m2 == n2 | m3 == n3 )
              H.rho = rho;
              pseudovar = pseudoinit(mol);
              %
              % compute the local ionic potential
              %
              [H.vion,vionT,rho1] = getvion(mol,pseudovar);
              %
              % compute nonlocal ionic potential
              %
              [H.wqmat,H.wqsign] = getwq(mol,pseudovar);
              clear pseudovar;
              % 
              % Calculate Hartree and exchange-correlation potential
              % from the known density rho
              % 
              [vhart,vxc,uxc2,rho]=getvhxc(mol,abs(rho));
              %
              % Save a copy of the Hartree and exchange-correlation
              % separately for DCM updating purposes
              %
              H.vnp = vhart+vxc;
              %
              % there may be other external potential
              % 
              vext = get(mol,'vext');
              H.vext = vext;
              %
              % vtot does not include non-local ionici potential
              % 
              H.vtot = getvtot(mol, H.vion, H.vext, vhart, vxc);
              H.rho  = rho;
              H.mpbc = get(mol,'mpbc');
              H.mpbc_denom = get(mol,'mpbc_denom');
              H.supercell = get(mol,'supercell');
              H = class(H,'Ham');
           else
              error('dimension of rho does not match that of the molecule');
           end; 
        else
           error('The 2nd input argument must be charge density (3D)');
        end;
     else
        error('The first input argument must be a Molecule object');
     end;
  otherwise
     error('Cannot have more than one arguement');
end;
