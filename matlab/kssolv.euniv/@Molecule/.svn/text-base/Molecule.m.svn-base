function mol = Molecule(varargin)
%
% purpose: molecule/crystal constructor
%
% usage: 1. one can use
%                 mol = Molecule();
%           to first define an empty Molecule object, and then
%           define various attributes of the Molecule (Bravais lattice, 
%           atoms, coordinates, energy cutoff etc.) by using the set method
%       
%        2. one can use
%                mol = Molecule(C);
%           to first specify the supercell C (a 3x3 matrix) of 
%           the molecule/crystal and then use the set method to define 
%           other attributes
%
%        3. one can use 
%                mol = Molecule(atomlist, xyzlist);
%           to first specify the list of atoms and their coordinates
%           and then use the set method to define other attributes
%
switch (nargin)
   case 0
     % set up a null molecule
     mol.name      = []; 
     mol.supercell = []; 
     mol.atomlist  = [];
     mol.atypes    = [];
     mol.xyzlist   = [];
     mol.vol       = [];
     mol.natoms    = [];
     mol.ecut      = [];
     mol.ecut2     = []; 
     mol.n1        = [];
     mol.n2        = [];
     mol.n3        = [];
     mol.rcut      = 3.8;
     mol.nel       = [];
     mol.nspin     = 1;
     mol.vext      = [];
     mol.temperature = 0.0;
     mol.mpbc      = 0;
     mol.mpbc_denom= 1;
     mol = class(mol,'Molecule');
   case 1
     %
     % the single input argument is Bravais lattice or energy cutoff or 
     % the name of the molecule
     %
     mol.name      = []; 
     mol.supercell = []; 
     mol.atomlist  = [];
     mol.atypes    = [];
     mol.xyzlist   = [];
     mol.vol       = [];
     mol.natoms    = [];
     mol.ecut      = [];
     mol.ecut2     = []; 
     mol.n1        = [];
     mol.n2        = [];
     mol.n3        = [];
     mol.nel       = [];
     mol.nspin     = 1;
     mol.rcut      = 3.8;
     mol.vext      = [];
     mol.temperature = 0.0;
     mol.mpbc      = 0;
     mol.mpbc_denom= 1;
     %
     % initialize the primitive (super)cell or name
     %
     if ( isnumeric(varargin{1}))
        [nrow, ncol] = size(varargin{1});
        if (nrow == 3 & ncol == 3)
           mol.supercell = varargin{1};
           mol.vol      = abs(det(varargin{1}));
           mol = class(mol,'Molecule');
        else
           error('The input arguement must be a 3x3 floating point array!');
        end;
     elseif ( ischar(varargin{1}) )
        mol.name = varargin{1};
        mol = class(mol,'Molecule');
     else
        error('The input arguement must be a 3x3 floating point array!');
     end;
   case 2
     %
     % the two input arguements must be atomlist and xyzlist 
     % (check the input type!)
     %
     atomlist = varargin{1};
     xyzlist  = varargin{2};
     if ( length(atomlist) == size(xyzlist,1) )
        mol.atomlist = atomlist;
        mol.atypes   = getatypes(atomlist);
        mol.xyzlist  = xyzlist;
        % set the Bravais lattice to the unit cube
        mol.name     = [];
        mol.supercell= [];
        mol.vol      = [];
        mol.natoms   = length(atomlist);
        mol.ecut     = [];
        mol.ecut2    = [];
        mol.n1       = [];
        mol.n2       = [];
        mol.n3       = [];
        mol.rcut     = 3.8;
        mol.nel      = getnel(mol);
        mol.nspin    = 1;
        mol.vext     = [];
        mol.temperature = 0.0;
        mol.mpbc      = 0;
        mol.mpbc_denom= 1;
        mol = class(mol,'Molecule');
     else
        error('The length of the atomlist must match that of the xyzlist');
     end
   case 3
     %
     % the three input arguements are: Bravais lattice, atomlist and xyzlist
     %
     mol.name     = [];
     mol.supercell= varargin{1};
     mol.atomlist = varargin{2};
     mol.atypes   = getatypes(varargin{2});
     mol.xyzlist  = varargin{3};
     mol.vol      = abs(det(varargin{1}));
     mol.natoms   = length(atomlist);
     mol.ecut     = []; 
     mol.ecut2    = [];
     mol.n1       = [];
     mol.n2       = [];
     mol.n3       = [];;
     mol.rcut     = 3.8;
     mol.nel      = getnel(mol);
     mol.nspin    = 1;
     mol.vext     = [];
     mol.temperature = 0.0;
     mol.mpbc      = 0;
     mol.mpbc_denom= 1;
     mol = class(mol,'Molecule');
   case 4
     %
     % the four input arguements are: Bravais lattice, atomlist and xyzlist
     % and the energy cutoff
     %
     mol.name     = [];
     C            = varargin{1};
     mol.supercell= C;
     mol.atomlist = varargin{2};
     mol.atypes   = getatypes(varargin{2});
     mol.xyzlist  = varargin{3};
     ecut         = varargin{4}; 
     mol.ecut     = ecut;
     mol.vol      = abs(det(C));
     mol.natoms   = length(atomlist);
     mol.ecut2    = 4*ecut;
     % we can estimate the minimum n1,n2,n3 values (based on sampling theorem)
     fackpt = 2*sqrt(2*ecut)/pi;
     d1 = fackpt*norm(C(:,1));
     d2 = fackpt*norm(C(:,2)); 
     d3 = fackpt*norm(C(:,3)); 
     mol.n1       = ceil(d1);
     mol.n2       = ceil(d2);
     mol.n3       = ceil(d3);
     mol.rcut     = 3.8;
     mol.nel      = getnel(mol);
     mol.nspin    = 1;
     mol.vext     = zeros(mol.n1,mol.n2,mol.n3);
     mol.temperature = 0.0;
     mol.mpbc      = 0;
     mol.mpbc_denom= 1;
     mol = class(mol,'Molecule');
   otherwise
     error('invalid Molecule constructor'); 
end;
