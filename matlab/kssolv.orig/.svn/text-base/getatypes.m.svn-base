function atypes = getatypes(atomlist)
%
%  Purpose: determine the number of different types of
%           atoms in the atomlist
%
%  Usage: atypes = getatypes(atomlist)l 
%
%  Input:
%     atomlist --- an array of Atoms objects that many contain repetitions. 
%
%  Output:
%     atypes   --- an array of Atoms objects that contain unique atoms
%
natoms = length(atomlist);
amask = zeros(110,1);
for i = 1:natoms
   anum = get(atomlist(i),'anum');
   amask(anum)=1;
end;
inz = find(amask~=0);
for j = 1:length(inz)
   atypes(j) = Atom(inz(j));
end;
