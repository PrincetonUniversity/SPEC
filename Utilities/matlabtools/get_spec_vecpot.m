function Acov = get_spec_vecpot(data,lvol,sarr,tarr,zarr)
 
%
% GET_SPEC_VECPOT( DATA, LVOL, SARR, TARR, ZARR )
% ===============================================
%
% Computes covariant components of A in volume lvol
%
% INPUT
% -----
%   -data     : must be produced by calling read_spec(filename)
%   -lvol     : is the volume number 
%   -sarr     : is the array of values for the s-coordinate
%   -tarr     : is the array of values for the theta-coordinate
%   -zarr     : is the array of values for the zeta-coordinate
%
% OUTPUT
% ------
%   -Acov     : cell structure with 2 arrays: A_theta, A_zeta each with size length(sarr)*length(tarr)*length(zarr)
%
% written by J.Loizu (2018)


Ate     = data.vector_potential.Ate{lvol};
Aze     = data.vector_potential.Aze{lvol};
Ato     = data.vector_potential.Ato{lvol};
Azo     = data.vector_potential.Azo{lvol};

Lrad    = data.input.physics.Lrad(lvol);

sarr    = transpose(sarr);
ns      = length(sarr);
nt      = length(tarr);
nz      = length(zarr);
sbar    = (sarr+1)/2;

mn      = data.output.mn;
im      = double(data.output.im);
in      = double(data.output.in);

At      = zeros(ns,nt,nz); % allocate data for vector potential along theta
Az      = zeros(ns,nt,nz); % allocate data for vector potential along zeta


% Construct Chebyshev polynomials 

T = get_spec_polynomial_basis(data, lvol, sarr);

% Construct regularization factors

fac = get_spec_regularization_factor(data, lvol, sarr, 'F');

% Construct vector potential covariant components

Lsingularity = false;
if (lvol==1) && (data.input.physics.Igeometry~=1)
  Lsingularity = true;
end


for l=1:Lrad+1
  for j=1:mn
    if( Lsingularity )
       basis = T{l}{1}(im(j));
      dbasis = T{l}{2}(im(j));
    else
       basis = T{l}{1};
      dbasis = T{l}{2};
    end

    for it=1:nt
      for iz=1:nz
       cosa = cos(im(j)*tarr(it)-in(j)*zarr(iz));
       sina = sin(im(j)*tarr(it)-in(j)*zarr(iz));
       At(:,it,iz) = At(:,it,iz) + fac{j}{1}.* basis.*( Ate(l,j)*cosa + Ato(l,j)*sina );
       Az(:,it,iz) = Az(:,it,iz) + fac{j}{1}.* basis.*( Aze(l,j)*cosa + Azo(l,j)*sina );
      end
    end
  end
end


Acov{1} = At;
Acov{2} = Az; 
