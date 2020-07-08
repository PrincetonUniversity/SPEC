function [sarr, out] = get_spec_R_derivatives(data, vol, sarr, tarr, zarr, RorZ)
%
% Returns the derivatives of R corresponding to SPEC coordinate 
% system. Stellarator symmetry assumed.
%
% INPUT
% -----
%	filename: 	SPEC hdf5 output filename
%	vol:		Volume number
%	sarr:       s-coordinate array, shape (ns, 1)
%	tarr:		Theta angle array
%	zarr:		Phi angle array
%   RorZ:       Returns either R derivatives (='R') or Z derivatives (='Z')
%
% OUTPUT
% ------
%	sarr:		s-coordinate array
%	R:		4xlength(sarr) array containing R, dR / ds, 
%			dR / dtheta and dR / dphi
%
% Written by A.Baillod (2019)
%


% Load geometry
mn     = data.mn;
im     = double(data.im);
in     = double(data.in);
Rmn = data.Rbc;
Zmn = data.Zbs;
    
% Allocate data for R and its derivative in s, theta and phi (4), for each
% and for ns points 
ns = length(sarr);
nt = length(tarr);
nz = length(zarr);
Rarr = cell(1,4); 
Zarr = cell(1,4);

for ii=1:4
   Rarr{ii} =  zeros(ns, nt, nz);
   Zarr{ii} =  zeros(ns, nt, nz);
end

% Compute the regularisation factor
factor = get_spec_regularisation_factor(data, vol, sarr, 'G');
                        
% And R derivatives
if RorZ=='R'
  for imn=1:mn
    for it=1:nt
      for iz=1:nz
        cosa = cos(double(im(imn)*tarr(it) - in(imn)*zarr(iz)));
        sina = sin(double(im(imn)*tarr(it) - in(imn)*zarr(iz)));

        Rarr{1}(:,it,iz) = Rarr{1}(:,it,iz) + ( Rmn(imn, vol  ) + (Rmn(imn, vol+1) - Rmn(imn, vol)) * factor{imn}{1})                   * cosa;
        Rarr{2}(:,it,iz) = Rarr{2}(:,it,iz) + (                   (Rmn(imn, vol+1) - Rmn(imn, vol)) * factor{imn}{2})                   * cosa;
        Rarr{3}(:,it,iz) = Rarr{3}(:,it,iz) - ( Rmn(imn, vol  ) + (Rmn(imn, vol+1) - Rmn(imn, vol)) * factor{imn}{1}) * double(im(imn)) * sina;
        Rarr{4}(:,it,iz) = Rarr{4}(:,it,iz) + ( Rmn(imn, vol  ) + (Rmn(imn, vol+1) - Rmn(imn, vol)) * factor{imn}{1}) * double(in(imn)) * sina;
      end
    end
  end
  out = Rarr;
  
elseif RorZ=='Z'
  for imn=1:mn
    for it=1:nt
      for iz=1:nz
        cosa = cos(double(im(imn)*tarr(it) - in(imn)*zarr(iz)));
        sina = sin(double(im(imn)*tarr(it) - in(imn)*zarr(iz)));

        Zarr{1}(:,it,iz) = Zarr{1}(:,it,iz) + ( Zmn(imn, vol  ) + (Zmn(imn, vol+1) - Zmn(imn, vol)) * factor{imn}{1})                  * sina;
        Zarr{2}(:,it,iz) = Zarr{2}(:,it,iz) + (                   (Zmn(imn, vol+1) - Zmn(imn, vol)) * factor{imn}{2})                  * sina;
        Zarr{3}(:,it,iz) = Zarr{3}(:,it,iz) + ( Zmn(imn, vol  ) + (Zmn(imn, vol+1) - Zmn(imn, vol)) * factor{imn}{1}) * double(im(imn))* cosa;
        Zarr{4}(:,it,iz) = Zarr{4}(:,it,iz) - ( Zmn(imn, vol  ) + (Zmn(imn, vol+1) - Zmn(imn, vol)) * factor{imn}{1}) * double(in(imn))* cosa;
      end
    end
  end    
  out = Zarr;
    
else
    error(['Not defined for RorZ=',RorZ]);
end
    

