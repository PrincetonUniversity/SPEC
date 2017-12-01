function plot_spec_grid(data,nz0,newfig)

% Produces plot of coordinate surfaces
%
% INPUT
%   -data     : must be produced by calling read_spec_grid(filename)
%   -nz0      : toroidal plane number at which coordinates are shown (nz0=1 at toroidal angle phi=0)
%   -newfig   : opens(=1) or not(=0) a new figure
%
%   written by J.Loizu (2015)


if(newfig==1)
figure
end

nvol   = data.Nvol;

Lrad   = data.Lrad;
Nt     = data.Nt;
Nz     = data.Nz;

Rij    = data.Rij;
Zij    = data.Zij;

ccol   = 'm';
cthick = 12;

iz     = nz0-1;

for i=1:nvol
  for l=1:Lrad(i)+1

    Rc   =  Rij(i,1+Nt*iz:(iz+1)*Nt,l);
    Zc   =  Zij(i,1+Nt*iz:(iz+1)*Nt,l);

    scatter(Rc,Zc,cthick,'filled',ccol)
    axis equal
    hold on
  end
end



