function Bcontrav = get_spec_magfield(fdata,lvol,sarr,tarr,zarr)
 
 
% Computes contravariant components of B in volume lvol
%
% INPUT
%   -fdata    : must be produced by calling read_spec_field(filename)
%   -lvol     : is the volume number 
%   -sarr     : is the array of values for the s-coordinate
%   -tarr     : is the array of values for the theta-coordinate
%   -zarr     : is the array of values for the zeta-coordinate
%
% OUTPUT
%   -Bcontrav : cell structure with 3 arrays: B^s, B^theta, B^zeta each with size length(sarr)*length(tarr)*length(zarr)
%
% Note: Stellarator symmetry is only assumed in the Jacobian evaluation
%
% written by J.Loizu (2016)


jac     = get_spec_jacobian(fdata,lvol,sarr,tarr,zarr);  % get jacobian of the coordinates (stell sym)

Ate     = fdata.Ate{lvol};
Aze     = fdata.Aze{lvol};
Ato     = fdata.Ato{lvol};
Azo     = fdata.Azo{lvol};

Lrad    = fdata.Lrad(lvol);

ns      = length(sarr);
nt      = length(tarr);
nz      = length(zarr);

mn      = fdata.mn;
im      = double(fdata.im);
in      = double(fdata.in);

Bs      = zeros(ns,nt,nz); % allocate data for magnetic field along s
Bt      = zeros(ns,nt,nz); % allocate data for magnetic field along theta
Bz      = zeros(ns,nt,nz); % allocate data for magnetic field along zeta

T       = cell(Lrad+1,2);  % allocate data for Chebyshev polynomials and their derivatives

T{1}{1} = ones(ns,1);
T{1}{2} = zeros(ns,1);

T{2}{1} = transpose(sarr);
T{2}{2} = ones(ns,1);


% Construct Chebyshev polynomials and their derivatives
for l=3:Lrad+1
  T{l}{1} = 2*transpose(sarr).*T{l-1}{1} - T{l-2}{1};
  T{l}{2} = 2*T{l-1}{1} + 2*transpose(sarr).*T{l-1}{2} - T{l-2}{2};
end

% Construct regularization factors and their derivatives
fac = get_spec_regularisation_factor(fdata, lvol, sarr, 'F');

% Construct magnetic field contravariant components
for l=1:Lrad+1
    for j=1:mn
        for it=1:nt
            for iz=1:nz
                cosa = cos(im(j)*tarr(it)-in(j)*zarr(iz));
                sina = sin(im(j)*tarr(it)-in(j)*zarr(iz));
                Bs(:,it,iz) = Bs(:,it,iz) +  fac{j}{1}.*T{l}{1}.*( (im(j)*Azo(l,j) + in(j)*Ato(l,j))*cosa - (im(j)*Aze(l,j) + in(j)*Ate(l,j))*sina  );
                Bt(:,it,iz) = Bt(:,it,iz) - (fac{j}{1}.*T{l}{2}+fac{j}{2}.*T{l}{1}).*( Aze(l,j)*cosa + Azo(l,j)*sina );
                Bz(:,it,iz) = Bz(:,it,iz) + (fac{j}{1}.*T{l}{2}+fac{j}{2}.*T{l}{1}).*( Ate(l,j)*cosa + Ato(l,j)*sina );
            end
        end
    end
end



Bcontrav{1} = Bs./jac; %actual B^s
Bcontrav{2} = Bt./jac; %actual B^theta
Bcontrav{3} = Bz./jac; %actual B^zeta
