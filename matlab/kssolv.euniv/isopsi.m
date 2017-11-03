function [h,psir] = isopsi(X,varargin)
%
%  usage: [h,psir] = isopsi(X,j);
%         [h,psir] = isopsi(X);
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
if (length(options) > 1)
   j = options{2};
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
   h = isosurface(psir);
   %daspect([1 1 1]); axis tight; 
   %colormap(cool)
   %camup([1 0 0 ]); campos([25 -55 5]) 
   %camlight; lighting phong
end;
h = gca;
