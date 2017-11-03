function display(atom)
%
% display method for atom
%
fprintf('atomic symbol: %s\n', atom.symbol);
fprintf('atomic number: %d\n', atom.anum);
fprintf('number of valence electrons: %d\n', atom.venum);
