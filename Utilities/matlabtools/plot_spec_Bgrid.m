function plot_spec_Bgrid(data,nz0,plotstyle,newfig)

% Produces plot of BR, Bphi, BZ, on the coordinate grid points
%
% INPUT
%   -data      : must be produced by calling read_spec_grid(filename)
%   -nz0       : toroidal plane number at which B is obtained (nz0=1 at toroidal angle phi=0)
%   -plotstyle : whether pcolor ('pcolor') or scatter ('scatter') is to be used
%   -newfig    : opens(=1) or not(=0) a new figure
%
% Note: The B components are the cylindrical contravariant terms, that is why Bphi is multiplied by R in order to get B \dot unitvectorphi.
%
% written by J.Loizu (2018)


if(newfig==1)
figure
end

nvol   = data.Mvol;

Lrad   = data.Lrad;
Nt     = data.Nt;
Nz     = data.Nz;
Nfp    = data.Nfp;

Rij    = data.Rij;
Zij    = data.Zij;
BR     = data.BR;  
Bp     = data.Bp; 
BZ     = data.BZ;  

iz     = nz0-1;
phi0   = double((2*pi/Nfp)*(iz/Nt));
rzdata = get_spec_rzarr(data,nvol,1,linspace(0,2*pi,32),phi0);
rmax   = max(rzdata{1});
rmin   = min(rzdata{1});
zmax   = max(rzdata{2});
zmin   = min(rzdata{2});

if(strcmp(plotstyle,'pcolor')==1)

for i=1:nvol
  ngrid  = Lrad(i)+1;

  if(i==1)
   nstart = 2;
  else
   nstart = 1;
  end

  np1  = Nt;
  np2  = 1+ngrid-nstart;

  Rc   = reshape(Rij(i,1+Nt*iz:(iz+1)*Nt,nstart:ngrid),np1,np2);
  Zc   = reshape(Zij(i,1+Nt*iz:(iz+1)*Nt,nstart:ngrid),np1,np2);
  br   = reshape(BR(i,1+Nt*iz:(iz+1)*Nt,nstart:ngrid),np1,np2);
  bp   = reshape(Bp(i,1+Nt*iz:(iz+1)*Nt,nstart:ngrid).*Rij(i,1+Nt*iz:(iz+1)*Nt,nstart:ngrid),np1,np2);
  bz   = reshape(BZ(i,1+Nt*iz:(iz+1)*Nt,nstart:ngrid),np1,np2);

  subplot(3,1,1)
  pcolor(Rc,Zc,br)
  axis equal; colorbar; hold on
  xlim([0.95*rmin 1.05*rmax])
  ylim([1.05*zmin 1.05*zmax])
  title('B_R')
  subplot(3,1,2)
  pcolor(Rc,Zc,bp)
  axis equal; colorbar; hold on
  xlim([0.95*rmin 1.05*rmax])
  ylim([1.05*zmin 1.05*zmax])
  title('B_{\phi}')
  subplot(3,1,3)
  pcolor(Rc,Zc,bz)
  axis equal; colorbar; hold on
  xlim([0.95*rmin 1.05*rmax])
  ylim([1.05*zmin 1.05*zmax])
  title('B_Z')
end

elseif(strcmp(plotstyle,'scatter')==1)

cthick = 12 ;

for i=1:nvol
  ngrid  = Lrad(i)+1;
  if(i==1)
   nstart = 2;
  else
   nstart = 1;
  end
  for l=nstart:ngrid
    Rc   =  Rij(i,1+Nt*iz:(iz+1)*Nt,l);
    Zc   =  Zij(i,1+Nt*iz:(iz+1)*Nt,l);
    colR =  BR(i,1+Nt*iz:(iz+1)*Nt,l);
    colp =  Bp(i,1+Nt*iz:(iz+1)*Nt,l).*Rc;
    colZ =  BZ(i,1+Nt*iz:(iz+1)*Nt,l);
    
    subplot(3,1,1)
    scatter(Rc,Zc,cthick,colR)
    axis equal; colorbar; hold on
    xlim([0.95*rmin 1.05*rmax])
    ylim([1.05*zmin 1.05*zmax])
    title('B_R')
    subplot(3,1,2)
    scatter(Rc,Zc,cthick,colp)
    axis equal; colorbar; hold on
    xlim([0.95*rmin 1.05*rmax])
    ylim([1.05*zmin 1.05*zmax])
    title('B_{\phi}')
    subplot(3,1,3)
    scatter(Rc,Zc,cthick,colZ)
    axis equal; colorbar; hold on
    xlim([0.95*rmin 1.05*rmax])
    ylim([1.05*zmin 1.05*zmax])
    title('B_Z')
  end
end


end

