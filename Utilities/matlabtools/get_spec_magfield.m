function Bcontrav = get_spec_magfield(data,lvol,sarr,tarr,zarr)
 
%
% GET_SPEC_MAGFIELD( DATA, LVOL, SARR, TARR, ZARR )
% =================================================
%
% Computes contravariant components of B in volume lvol
%
% INPUT
% -----
%   -data    : must be produced by calling read_spec(filename)
%   -lvol     : is the volume number 
%   -sarr     : is the array of values for the s-coordinate
%   -tarr     : is the array of values for the theta-coordinate
%   -zarr     : is the array of values for the zeta-coordinate
%
% OUTPUT
% ------
%   -Bcontrav : cell structure with 3 arrays: B^s, B^theta, B^zeta each with size length(sarr)*length(tarr)*length(zarr)
%
% Note: Stellarator symmetry is only assumed in the Jacobian evaluation
%
% written by J.Loizu (2016)


jac     = get_spec_jacobian(data,lvol,sarr,tarr,zarr);  % get jacobian of the coordinates (stell sym)

Ate     = data.vector_potential.Ate{lvol};
Aze     = data.vector_potential.Aze{lvol};
Ato     = data.vector_potential.Ato{lvol};
Azo     = data.vector_potential.Azo{lvol};

Lrad    = data.input.physics.Lrad(lvol);

ns      = length(sarr);
nt      = length(tarr);
nz      = length(zarr);

mn      = data.output.mn;
im      = double(data.output.im);
in      = double(data.output.in);

Bs      = zeros(ns,nt,nz); % allocate data for magnetic field along s
Bt      = zeros(ns,nt,nz); % allocate data for magnetic field along theta
Bz      = zeros(ns,nt,nz); % allocate data for magnetic field along zeta

% Construct Chebyshev polynomials and their derivatives

T = get_spec_polynomial_basis(data, lvol, sarr');

% Construct magnetic field contravariant components
Lsingularity = false;
if (lvol==1) && (data.input.physics.Igeometry~=1)
  Lsingularity = true;
end

for l=1:Lrad+1
    for j=1:mn
        if Lsingularity
          basis  = T{l}{1}(im(j)+1,:);
          dbasis = T{l}{2}(im(j)+1,:); 
          basis = basis';
          dbasis = dbasis';
        else
	      basis  = T{l}{1};
          dbasis = T{l}{2};
	    end

        for it=1:nt
            for iz=1:nz
                cosa = cos(im(j)*tarr(it)-in(j)*zarr(iz));
                sina = sin(im(j)*tarr(it)-in(j)*zarr(iz));
                Bs(:,it,iz) = Bs(:,it,iz) +  basis .* ( (im(j)*Azo(l,j) + in(j)*Ato(l,j))*cosa - (im(j)*Aze(l,j) + in(j)*Ate(l,j))*sina  );
                Bt(:,it,iz) = Bt(:,it,iz) - dbasis .* ( Aze(l,j)*cosa + Azo(l,j)*sina );
                Bz(:,it,iz) = Bz(:,it,iz) + dbasis .* ( Ate(l,j)*cosa + Ato(l,j)*sina );
            end
        end
    end
end



Bcontrav{1} = Bs./jac; %actual B^s
Bcontrav{2} = Bt./jac; %actual B^theta
Bcontrav{3} = Bz./jac; %actual B^zeta
