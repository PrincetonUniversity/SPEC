function plot_spec_Bgrid(data,nz0,plotstyle,newfig)

%
% PLOT_SPEC_BGRID( DATA, NZ0, PLOTSTYLE, NEWFIG )
% ===============================================
%
% Produces plot of BR, Bphi, BZ, on the coordinate grid points
%
% INPUT
% -----
%   -data      : must be produced by calling read_spec(filename)
%   -nz0       : toroidal plane number at which B is obtained (nz0=1 at toroidal angle phi=0)
%   -plotstyle : whether pcolor ('pcolor') or scatter ('scatter') is to be used
%   -newfig    : opens(=1) or not(=0) a new figure
%
% Note: The B components are the cylindrical contravariant terms, that is why Bphi is multiplied by R in order to get B \dot unitvectorphi.
%
% written by J.Loizu (2018)
%
% OUTDATED - NEED DEBUG

% Check inputs
if ~strcmp(string(plotstyle),string('pcolor')) && ~strcmp(string(plotstyle),string('scatter'))
    error('InputError: Invalid plotstyle')
end
if nz0<1
    error('nzo should be greater than zero')
end
switch newfig
    case 0
        hold on
    case 1
        figure('Color','w','Position',[200 200 1500, 600])
        hold on
    case 2
        hold off
    otherwise
        error('Invalide newfig')
end

Mvol   = data.output.Mvol;

Lrad   = data.input.physics.Lrad;
Nt     = data.grid.Nt;
Nfp    = data.input.physics.Nfp;


iz     = nz0-1;
phi0   = double((2*pi/Nfp)*(iz/Nt));
R      = get_spec_R_derivatives(data, Mvol, 1, linspace(0,2*pi,64), phi0, 'R');
Z      = get_spec_R_derivatives(data, Mvol, 1, linspace(0,2*pi,64), phi0, 'Z');

rmax   = max(R{1});
rmin   = min(R{1});
zmax   = max(Z{1});
zmin   = min(Z{1});

if(strcmp(plotstyle,'pcolor')==1)

    for i=1:Mvol
      ngrid  = Lrad(i)+1;

      % Read data corresponding to correct volume
      Rij = data.grid.Rij{i};
      Zij = data.grid.Zij{i};
      BR  = data.grid.BR{i};  
      Bp  = data.grid.Bp{i}; 
      BZ  = data.grid.BZ{i};  

      % Reshape as an array
      if(i==1)
       nstart = 2;
      else
       nstart = 1;
      end

      np1  = Nt;
      np2  = 1+ngrid-nstart;

      Rc   = reshape(Rij(1+Nt*iz:(iz+1)*Nt,nstart:ngrid),np1,np2);
      Zc   = reshape(Zij(1+Nt*iz:(iz+1)*Nt,nstart:ngrid),np1,np2);
      br   = reshape(BR( 1+Nt*iz:(iz+1)*Nt,nstart:ngrid),np1,np2);
      bp   = reshape(Bp( 1+Nt*iz:(iz+1)*Nt,nstart:ngrid) ...
                   .*Rij(1+Nt*iz:(iz+1)*Nt,nstart:ngrid),np1,np2);
      bz   = reshape(BZ( 1+Nt*iz:(iz+1)*Nt,nstart:ngrid),np1,np2);

      % Double first entry to fill entire plane
      Rc(end+1,:) = Rc(1,:);
      Zc(end+1,:) = Zc(1,:);
      br(end+1,:) = br(1,:);
      bp(end+1,:) = bp(1,:);
      bz(end+1,:) = bz(1,:);

      % Plots
      subplot(1,3,1)
      pcolor(Rc,Zc,br)
      shading interp
      axis equal; colorbar; hold on
      xlim([0.95*rmin 1.05*rmax])
      ylim([1.05*zmin 1.05*zmax])
      title('B_R')

      subplot(1,3,2)
      pcolor(Rc,Zc,bp)
      shading interp
      axis equal; colorbar; hold on
      xlim([0.95*rmin 1.05*rmax])
      ylim([1.05*zmin 1.05*zmax])
      title('B_{\phi}')

      subplot(1,3,3)
      pcolor(Rc,Zc,bz)
      shading interp
      axis equal; colorbar; hold on
      xlim([0.95*rmin 1.05*rmax])
      ylim([1.05*zmin 1.05*zmax])
      title('B_Z')
    end

elseif(strcmp(plotstyle,'scatter')==1)

    cthick = 12 ;

    for i=1:Mvol
      % Read data corresponding to correct volume
      Rij = data.grid.Rij{i};
      Zij = data.grid.Zij{i};
      BR  = data.grid.BR{i};  
      Bp  = data.grid.Bp{i}; 
      BZ  = data.grid.BZ{i};  
      
      ngrid  = Lrad(i)+1;
      if(i==1)
       nstart = 2;
      else
       nstart = 1;
      end
      for l=nstart:ngrid
        Rc   =  Rij(1+Nt*iz:(iz+1)*Nt,l);
        Zc   =  Zij(1+Nt*iz:(iz+1)*Nt,l);
        colR =  BR(1+Nt*iz:(iz+1)*Nt,l);
        colp =  Bp(1+Nt*iz:(iz+1)*Nt,l).*Rc;
        colZ =  BZ(1+Nt*iz:(iz+1)*Nt,l);

        subplot(1,3,1)
        scatter(Rc,Zc,cthick,colR)
        axis equal; colorbar; hold on
        xlim([0.95*rmin 1.05*rmax])
        ylim([1.05*zmin 1.05*zmax])
        title('B_R')
        
        subplot(1,3,2)
        scatter(Rc,Zc,cthick,colp)
        axis equal; colorbar; hold on
        xlim([0.95*rmin 1.05*rmax])
        ylim([1.05*zmin 1.05*zmax])
        title('B_{\phi}')
        
        subplot(1,3,3)
        scatter(Rc,Zc,cthick,colZ)
        axis equal; colorbar; hold on
        xlim([0.95*rmin 1.05*rmax])
        ylim([1.05*zmin 1.05*zmax])
        title('B_Z')
      end
    end


end

