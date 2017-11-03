function X0 = genX0(mol,ncols)
%
% Usage: X0 = genX0(mol,ncols); 
%        X0 = genX0(mol); 
%
% Purpose:
%    Generate the initial wavefunction
%
% Input:
%    mol --- a Molecule object
%    ncols --- (optional) the number of wavefunctions included
% Output:
%    X0  --- a Wavefun object. The object stores the Fourier
%            coefficients associated with freqency (wave) vectors
%            that satifies the kinetic energy cutoff defined in 
%            the mol object.  These frequency vectors 
%            are generated and stored in FreqMask(mol).  
%
nspin = get(mol,'nspin');
if (nargin == 1)
   ncols = getnel(mol)/2*nspin;
end;

n1   = get(mol,'n1');
n2   = get(mol,'n2');
n3   = get(mol,'n3');
n123 = n1*n2*n3;

gmask = FreqMask(mol);
gm    = get(gmask,'gm');
idxnz = find(gm~=0);

W = [];   
for j = 1:ncols
   psir = randn(n1,n2,n3);
   psif = fftn(psir);
   W = [W psif(idxnz)];
end;
[Q,R]=qr(W,0);
X0 = Wavefun(Q,n1,n2,n3,idxnz);
