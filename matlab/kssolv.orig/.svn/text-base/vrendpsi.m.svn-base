function [h,psir] = vrendpsi(X,varargin)
%
%  usage: [h,psir] = vrendpsi(X,j);
%         [h,psir] = vrendpsi(X);
%
%  view the isosurface of of a wavefunction
%
%  input:
%    X    -  a Wavefun that contains one or more 3-D wavefunctions
%    j    -  index of the wavefunction to be displayed. If no index
%            is used, the first wavefunction will be displayed.
%  output:
%    h    -  graphic handle of the displayed wavefunction
%    psir - wavefunction stored in a 3-D array
%
clf;
options = varargin;
nopt = length(options)
if ( nopt > 0)
   j = options{1};
else 
   j = 1;
end;
%
psi = get(X,'psi');
%
% find out the number of wavefunctions
%
numpsi = size(psi,2);
if (j < 1 || j > numpsi)
   fprintf('The second argument must be less than %d\n',numpsi);
   return;
else
   psir = abs(fftshift(ifftn(psi{j})));
   h = vol3d('cdata',psir,'texture','2D');
   view(3); 
   % Update view since 'texture' = '2D'
   vol3d(h);  
   alphamap('rampdown'), alphamap('decrease'), alphamap('decrease')
   alphamap('rampup');
   %alphamap(.06 .* alphamap);
end;
