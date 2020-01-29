function plot_spec_modB_boundary(data,Nvol,nt,nz)

% Produces plot of |B| on the full boundary surface in toroidal geometry.
%
% INPUT
%   -data   : data obtained via read_spec(filename)
%   -Nvol   : total number of volumes
%   -nt     : poloidal resolution for the plotting (e.g. nt=64)
%   -nz     : toroidal resolution for the plotting (e.g. nz=64)
%
% written by J.Loizu (2016)
% modified by J.Loizu (01.2017)

sarr = 1;

tarr = linspace(0,2*pi,nt);

zarr = linspace(0,2*pi,nz);


% Read vector potential

fdata  = fdata_from_data(data);


% Compute |B|

modB   = get_spec_modB(fdata,Nvol,sarr,tarr,zarr);


% Compute function (R,Z)(s,theta,zeta)

rzdata = get_spec_rzarr(fdata,Nvol,sarr,tarr,zarr);

R = squeeze(rzdata{1});   
Z = squeeze(rzdata{2});


% Construct cartesian corrdinates 

X = zeros(nt,nz);

Y = zeros(nt,nz);
 
for it=1:nt
 for iz=1:nz
  X(it,iz) = R(it,iz)*cos(zarr(iz));
  Y(it,iz) = R(it,iz)*sin(zarr(iz));
 end
end
 
 
% Plot

figure

h=surf(X,Y,Z,squeeze(modB(1,:,:)));

axis equal
shading interp
colorbar
title('| B |')
