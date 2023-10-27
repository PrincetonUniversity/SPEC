function plot_spec_grid(data,nz0,newfig)

% 
% PLOT_SPEC_GRID( DATA, NZ0, NEWFIG )
% ===================================
%
% Produces plot of coordinate surfaces
%
% INPUT
% -----
%   -data     : must be produced by calling read_spec(filename)
%   -nz0      : toroidal plane number at which coordinates are shown (nz0=1 at toroidal angle phi=0)
%   -newfig   : opens(=1) or not(=0) a new figure, or overwrite selected
%               figure (=2)
%
%   written by J.Loizu (2015)
%

switch newfig
    case 0
        hold on
    case 1
        figure
        hold on
    case 2
        hold off
    otherwise
        error('InputError: invalid newfig')
end

Mvol   = data.output.Mvol;

Lrad   = data.input.physics.Lrad;
Nt     = data.grid.Nt;

Rij    = data.grid.Rij;
Zij    = data.grid.Zij;

ccol   = 'm';
cthick = 6;

iz     = nz0-1;

for i=1:Mvol
  for l=1:Lrad(i)+1

    R_tmp = Rij{i};
    Z_tmp = Zij{i};
      
    Rc   =  R_tmp(1+Nt*iz:(iz+1)*Nt,l);
    Zc   =  Z_tmp(1+Nt*iz:(iz+1)*Nt,l);

    scatter(Rc,Zc,cthick,'filled',ccol)
    axis equal
    hold on
  end
end



