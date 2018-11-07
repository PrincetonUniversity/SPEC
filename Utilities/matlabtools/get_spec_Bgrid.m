function bdata = get_spec_Bgrid(data,nz0,lvol)

% Obtains the canonical cylindrical components of B, namely (B^R, R*B^phi, B^Z), on the coordinate grid points
%
% INPUT
%   -data     : must be produced by calling read_spec_grid(filename)
%   -nz0      : toroidal plane number at which B is obtained (nz0=1 at toroidal angle phi=0)
%   -lvol     : volume number in which B is obtained
%
% OUTPUT
%   -bdata    : cell of size 5 containing the values of R, Z, B^R, R*B^phi, B^Z on the grid points
%
% written by J.Loizu (2018)

Lrad   = data.Lrad;
Nt     = data.Nt;
Nz     = data.Nz;

Rij    = data.Rij;
Zij    = data.Zij;
BR     = data.BR;  
Bp     = data.Bp; 
BZ     = data.BZ;  

iz     = nz0-1;
ngrid  = Lrad(lvol)+1;

if(lvol==1)
   nstart = 2;
else
   nstart = 1;
end

bdata{1} =  squeeze(Rij(lvol,1+Nt*iz:(iz+1)*Nt,nstart:ngrid));
bdata{2} =  squeeze(Zij(lvol,1+Nt*iz:(iz+1)*Nt,nstart:ngrid));
bdata{3} =  squeeze(BR(lvol,1+Nt*iz:(iz+1)*Nt,nstart:ngrid));
bdata{4} =  squeeze(Bp(lvol,1+Nt*iz:(iz+1)*Nt,nstart:ngrid)).*bdata{1};
bdata{5} =  squeeze(BZ(lvol,1+Nt*iz:(iz+1)*Nt,nstart:ngrid));


