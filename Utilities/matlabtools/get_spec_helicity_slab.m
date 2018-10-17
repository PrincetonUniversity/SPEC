function H = get_spec_helicity_slab(fdata,lvol,ns,nt,nz)
 
 
% Calculates Magnetic Helicity in a given volume in slab geometry
%
% INPUT
%   -data      : must be produced by calling read_spec_field(filename)
%   -lvol      : volume in which the helicity is evaluated
%   -ns        : is the s-coordinate resolution 
%   -nt        : is the theta-coordinate resolution
%   -nz        : is the zeta-coordinate resolution
%
% OUTPUT
%   -H         : Magnetic Helicity in volume lvol 
%
% Note: Stellarator symmetry is assumed
% Note: This definition of Helicity, integral(A.B), is not gauge-invariant for multi-valued gauges
%
% written by J.Loizu (2018)

Ate     = fdata.Ate{lvol};
Aze     = fdata.Aze{lvol};

sarr    = linspace(-1,1,ns);
tarr    = linspace(0,2*pi,nt);
zarr    = linspace(0,2*pi,nz);
sarr    = transpose(sarr);
sbar    = (sarr+1)/2;

Lrad    = fdata.Lrad(lvol);
mn      = fdata.mn;
im      = double(fdata.im);
in      = double(fdata.in);

h1      = zeros(ns,nt,nz); % allocate data for magnetic helicity integrand 1
h2      = zeros(ns,nt,nz); % allocate data for magnetic helicity integrand 2
h3      = zeros(ns,nt,nz); % allocate data for magnetic helicity integrand 3
h4      = zeros(ns,nt,nz); % allocate data for magnetic helicity integrand 4
T       = cell(Lrad+1,2);  % allocate data for Chebyshev polynomials and their derivatives
fac     = cell(mn,2);      % allocate data for regularization factors and their derivatives

T{1}{1} = ones(ns,1);
T{1}{2} = zeros(ns,1);

T{2}{1} = sarr;
T{2}{2} = ones(ns,1);


% Construct Chebyshev polynomials and their derivatives

for l=3:Lrad+1
  T{l}{1} = 2*sarr.*T{l-1}{1} - T{l-2}{1};
  T{l}{2} = 2*T{l-1}{1} + 2*sarr.*T{l-1}{2} - T{l-2}{2};
end


% Construct regularization factors and their derivatives

for j=1:mn
  fac{j}{1}  = ones(ns,1);
  fac{j}{2}  = zeros(ns,1);
end


% Construct magnetic helicity integrand

for l=1:Lrad+1
  for j=1:mn
    for it=1:nt
      for iz=1:nz
       cosa = cos(im(j)*tarr(it)-in(j)*zarr(iz));
       sina = sin(im(j)*tarr(it)-in(j)*zarr(iz));
       h1(:,it,iz) = h1(:,it,iz) + fac{j}{1}.*T{l}{1}.*Ate(l,j)*cosa; %A_t
       h2(:,it,iz) = h2(:,it,iz) - (fac{j}{1}.*T{l}{2}+fac{j}{2}.*T{l}{1}).*Aze(l,j)*cosa;  % -dA_z/ds
       h3(:,it,iz) = h3(:,it,iz) + fac{j}{1}.*T{l}{1}.*Aze(l,j)*cosa; %A_z
       h4(:,it,iz) = h4(:,it,iz) + (fac{j}{1}.*T{l}{2}+fac{j}{2}.*T{l}{1}).*Ate(l,j)*cosa;  % +dA_t/ds
      end
    end
  end
end

h = h1.*h2 + h3.*h4;

H = trapz(sarr,trapz(tarr,trapz(zarr,h,3),2));


