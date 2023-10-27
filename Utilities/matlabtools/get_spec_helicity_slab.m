function H = get_spec_helicity_slab(data,lvol,ns,nt,nz)
 
% 
% GET_SPEC_HELICITY_SLAB( DATA, LVOL, NS, NT, NZ )
% ================================================
%
% Calculates Magnetic Helicity in a given volume in slab geometry
%
% INPUT
% -----
%   -data      : must be produced by calling read_spec(filename)
%   -lvol      : volume in which the helicity is evaluated
%   -ns        : is the s-coordinate resolution 
%   -nt        : is the theta-coordinate resolution
%   -nz        : is the zeta-coordinate resolution
%
% OUTPUT
% ------
%   -H         : Magnetic Helicity in volume lvol 
%
% Note: Stellarator symmetry is assumed
% Note: This definition of Helicity, integral(A.B), is not gauge-invariant for multi-valued gauges
%
% written by J.Loizu (2018)

    % Test input
    Igeometry = data.input.physics.Igeometry;
    if( Igeometry~=1 )
        error('Invalid geometry. Only available for slab')
    end

    % Read data
    Ate     = data.vector_potential.Ate{lvol};
    Aze     = data.vector_potential.Aze{lvol};

    sarr    = linspace(-1,1,ns);
    tarr    = linspace(0,2*pi,nt);
    zarr    = linspace(0,2*pi,nz);
    sarr    = transpose(sarr);

    Lrad    = data.input.physics.Lrad(lvol);
    mn      = data.output.mn;
    im      = double(data.output.im);
    in      = double(data.output.in);

    h1      = zeros(ns,nt,nz); % allocate data for magnetic helicity integrand 1
    h2      = zeros(ns,nt,nz); % allocate data for magnetic helicity integrand 2
    h3      = zeros(ns,nt,nz); % allocate data for magnetic helicity integrand 3
    h4      = zeros(ns,nt,nz); % allocate data for magnetic helicity integrand 4

    % Construct Chebyshev polynomials and their derivatives

    T = get_spec_polynomial_basis(data, lvol, sarr);

    % Construct regularization factors and their derivatives

    fac = get_spec_regularisation_factor(data,lvol,sarr,'F');

    % Construct magnetic helicity integrand

    Lsingularity = false;
    if (lvol==1) && (data.input.physics.Igeometry~=1)
      Lsingularity = true;
    end

    for l=1:Lrad+1
      for j=1:mn
        if Lsingularity
          basis  = T{l}{1}(im(j)+1);
          dbasis = T{l}{2}(im(j)+1);
        else
          basis  = T{l}{1};
          dbasis = T{l}{2};
        end

        for it=1:nt
          for iz=1:nz
           cosa = cos(im(j)*tarr(it)-in(j)*zarr(iz));
           h1(:,it,iz) = h1(:,it,iz) +  fac{j}{1}.* basis                   .*Ate(l,j)*cosa; %A_t
           h2(:,it,iz) = h2(:,it,iz) - (fac{j}{1}.*dbasis+fac{j}{2}.* basis).*Aze(l,j)*cosa;  % -dA_z/ds
           h3(:,it,iz) = h3(:,it,iz) +  fac{j}{1}.* basis                   .*Aze(l,j)*cosa; %A_z
           h4(:,it,iz) = h4(:,it,iz) + (fac{j}{1}.*dbasis+fac{j}{2}.* basis).*Ate(l,j)*cosa;  % +dA_t/ds
          end
        end
      end
    end

    % Evaluate helicity
    h = h1.*h2 + h3.*h4;
    H = trapz(sarr,trapz(tarr,trapz(zarr,h,3),2));

end
