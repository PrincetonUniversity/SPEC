function plot_spec_outfluxfun(data,ns,nt,z0,ncont,newfig)

%
% PLOT_SPEC_OUTFLUXFUN( DATA, NS, NT, Z0, NCONT, NEWFIG )
% =======================================================
%
% Plots iso-contours of Az on a given cross-section
%
% INPUT
% -----
%   -data     : must be produced by calling read_spec(filename)
%   -ns       : radial resolution for construction of Az
%   -nt       : poloidal resolution for construction of Az
%   -z0       : toroidal angle at which Az is evaluated
%   -ncont    : number of iso-contour lines
%   -newfig   : flag for whether a new figure should be open (=1) or not(=0)
%
% Note: should only work for free-boundary equilibria
% Note: poloidal flux function is psi=-A_phi if A_phi is the covariant component and psi=-R*A_phi if A_phi is the canonical component
%
% written by J.Loizu (2018)


Nvol = data.input.physics.Nvol;

lvol = Nvol+1;

sarr = linspace(-1,1,ns);

tarr = linspace(0,2*pi,nt);

acov = get_spec_vecpot(data,lvol,sarr,tarr,z0);
rz   = get_spec_rzarr( data,lvol,sarr,tarr,z0);

ffun = -acov{2}; 

if(newfig==1)
 figure; hold on; 
end

contour(rz{1},rz{2},ffun,ncont,'k')
