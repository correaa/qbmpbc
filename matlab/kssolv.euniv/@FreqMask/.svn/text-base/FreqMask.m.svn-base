function Gmask = FreqMask(varargin)
%
% FreqMask (Frequency Mask) constructor
% usage:  Gmask = FreqMask(mol);
%         Gmask = FreqMask(mol,ecut);
%
% Keep desired frequency vectors in gkk, gkx, gky, gkz, all others
% are set to zeros.
%
% If the second argument (energy cutoff ecut) is not specified,
% FreqMask will use the ecut defined in mol as the cutoff
% frequency.
%
if (nargin == 0)
   Gmask.ecut= [];
   Gmask.gm  = [];
   Gmask.ng  = [];
   Gmask.gkk = [];
   Gmask.gkx = [];
   Gmask.gky = [];
   Gmask.gkz = [];
   Gmask  = class(Gmask, 'FreqMask');
elseif (nargin < 3)
   switch (nargin)
      case 1
         %
         % mol is the only input
         %
         if ( isa(varargin{1},'Molecule') )
            %
            % check the validity of the Molecule
            %
            % ....
            mol = varargin{1};  
            n1 = get(mol,'n1');
            n2 = get(mol,'n2');
            n3 = get(mol,'n3');
            C = get(mol,'supercell');
            ecut = get(mol,'ecut');
            Gmask.ecut = ecut;
            ecut = ecut/2; 
            mpbc = get(mol,'mpbc');
         else
            error('The input argument must be a Molecule');
         end;
      case 2
         if ( isa(varargin{1},'Molecule') & isnumeric(varargin{2}) )
            %
            % check the validity of the Molecule
            %
            % ....
            mol = varargin{1};  
            n1 = get(mol,'n1');
            n2 = get(mol,'n2');
            n3 = get(mol,'n3');
            C = get(mol,'supercell');
            ecut = varargin{2};
            Gmask.ecut = ecut;
            ecut = ecut/2; 
            mpbc = get(mol,'mpbc');
         else
            error('The input argument must be a Molecule');
         end;
   end; 

   Gmask.gm   = zeros(n1,n2,n3);
   Gmask.ng   = 0;
   Gmask.gkk  = zeros(n1,n2,n3);
   Gmask.gkx  = zeros(n1,n2,n3);
   Gmask.gky  = zeros(n1,n2,n3);
   Gmask.gkz  = zeros(n1,n2,n3);

   CI = inv(C);

   for k3=1:n3
      k1=k3-1;
      if (k1 > n3/2) 
         k1=k1-n3;
      end
      for j3=1:n2
         j1=j3-1;
         if(j1 > n2/2) 
            j1=j1-n2;
         end;
         for i3=1:n1
            i1=i3-1;
            if (i1 > n1/2) 
               i1=i1-n1;
            end;  
            akkx=2*pi*(CI(1,1)*i1+CI(1,2)*j1+CI(1,3)*k1);
            akky=2*pi*(CI(2,1)*i1+CI(2,2)*j1+CI(2,3)*k1);
            akkz=2*pi*(CI(3,1)*i1+CI(3,2)*j1+CI(3,3)*k1);

            akk=(akkx^2+akky^2+akkz^2)/2;

            if (mpbc == 0 && akk <= ecut) 
%               Gmask.gm(k3,j3,i3)  = 1;
%               Gmask.gkk(k3,j3,i3) = akk;
%               Gmask.gkx(k3,j3,i3) = akkx;
%               Gmask.gky(k3,j3,i3) = akky;
%               Gmask.gkz(k3,j3,i3) = akkz;
               Gmask.gm(i3,j3,k3)  = 1;
               Gmask.gkk(i3,j3,k3) = akk;
               Gmask.gkx(i3,j3,k3) = akkx;
               Gmask.gky(i3,j3,k3) = akky;
               Gmask.gkz(i3,j3,k3) = akkz;
               Gmask.ng=Gmask.ng+1;
            elseif (mpbc > 0)
               Gmask.gm(i3,j3,k3)  = 1;
               Gmask.gkk(i3,j3,k3) = akk;
               Gmask.gkx(i3,j3,k3) = akkx;
               Gmask.gky(i3,j3,k3) = akky;
               Gmask.gkz(i3,j3,k3) = akkz;
               Gmask.ng=Gmask.ng+1;                
            end;
         end;
      end;
   end;
   Gmask = class(Gmask,'FreqMask');
else
   error('The number of arguments must be less than 3');
end; 
