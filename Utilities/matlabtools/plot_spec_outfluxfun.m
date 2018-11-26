function plot_spec_outfluxfun(fdata,ns,nt,z0,ncont,newfig)

 
% Plots iso-contours of -R*Az on a given cross-section
%
% INPUT
%   -fdata    : must be produced by calling read_spec_field(filename)
%   -ns       : radial resolution for construction of Az
%   -nt       : poloidal resolution for construction of Az
%   -z0       : toroidal angle at which Az is evaluated
%   -ncont    : number of iso-contour lines
%   -newfig   : flag for whether a new figure should be open (=1) or not(=0)
%
% Note: should only work for free-boundary equilibria
%
% written by J.Loizu (2018)


Nvol = fdata.Nvol;

lvol = Nvol+1;

sarr = linspace(-1,1,ns);

tarr = linspace(0,2*pi,nt);

acov = get_spec_vecpot(fdata,lvol,sarr,tarr,z0);

rz   = get_spec_rzarr(fdata,lvol,sarr,tarr,z0);

ffun = -rz{1}.*acov{2}; % poloidal flux function is psi=-R*A_phi

if(newfig==1)
 figure; hold on; 
end

contour(rz{1},rz{2},ffun,ncont,'k')
