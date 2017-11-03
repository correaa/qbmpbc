function ierr = finalize(a)
%
% Purpose: perform some error checking for a Molecule object
%
% usage: ierr = finalize(a);
%        
% Input:  a (Molecule) --- A Molecule object
% Output: ierr (integer) --- error flag. ierr=0 no error
%                                        ierr=1 the volume of the molecule
%                                        is larger than the unit cell volume
%
ierr = 0;
C = a.supercell;
Cvol = abs(det(C));
delx = max(a.xyzlist(:,1))-min(a.xyzlist(:,1));
dely = max(a.xyzlist(:,2))-min(a.xyzlist(:,2));
delz = max(a.xyzlist(:,3))-min(a.xyzlist(:,3));
Mvol = delx*dely*delz;
if (Cvol < Mvol)
   ierr = 1;
end;
