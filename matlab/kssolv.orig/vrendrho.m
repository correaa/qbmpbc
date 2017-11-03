function h = vrendrho(rho)
%
%  usage: h = vrendrho(rho);
%
%  view the isosurface of a charge density
%
%  input:
%    rho  -  a 3-D array to be viewed.
%  output:
%    h    -  graphic handle of the displayed rho
%
clf;
h = vol3d('cdata',rho,'texture','2D');
view(3); 
% Update view since 'texture' = '2D'
vol3d(h);  
alphamap('rampdown'), alphamap('decrease'), alphamap('decrease')
alphamap('rampup');
%alphamap(.06 .* alphamap);
