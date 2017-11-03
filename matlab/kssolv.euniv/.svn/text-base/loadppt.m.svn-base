function vwr = loadppt(atom)
%
% Usage: vwr = loadppt(atom);
% 
% Purpose:
%    Load the pseudo potential information associated with the atom object
%
% Input:
%   atom   --- an Atom object
%   vwr    --- pseudopotential related quantities (matrix)

% check the input argument
%
symbol = get(atom,'symbol');
kssolvpath = getenv('KSSOLVPATH');
if isempty(kssolvpath)
    kssolvpath = '.';
end
pseudopotpath = [kssolvpath '/PseudoPot'];
name   = sprintf('%s/%s',pseudopotpath, symbol);
pseudo = load(name);
vwr = pseudo.vwr;
