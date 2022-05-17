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

    % Test input
    Mvol = data.output.Mvol;
    if lvol<1 || lvol>Mvol
        error('InputError: invalid lvol')
    end

    if isempty(sarr)
        error('InputError: empty sarr')
    end
    if isempty(tarr)
        error('InputError: empty tarr')
    end
    if isempty(zarr)
        error('InputError: empty zarr')
    end

    if sarr(1)<-1 || sarr(end)>1
        error('InputError: invalid sarr')
    end

    % Read data 
    Ate     = data.vector_potential.Ate{lvol};
    Aze     = data.vector_potential.Aze{lvol};
    Ato     = data.vector_potential.Ato{lvol};
    Azo     = data.vector_potential.Azo{lvol};

    Lrad    = data.input.physics.Lrad(lvol);

    if(size(sarr,1)==1)
        sarr    = transpose(sarr);
    end
    ns      = length(sarr);
    nt      = length(tarr);
    nz      = length(zarr);

    mn      = data.output.mn;
    im      = double(data.output.im);
    in      = double(data.output.in);

    At      = zeros(ns,nt,nz); % allocate data for vector potential along theta
    Az      = zeros(ns,nt,nz); % allocate data for vector potential along zeta

    % Construct vector potential covariant components
    T = get_spec_polynomial_basis(data, lvol, sarr);
    
    % Construct vector potential covariant components
    Lsingularity = (lvol==1) && (data.input.physics.Igeometry~=1);

    for l=1:Lrad+1
        for j=1:length(im)
            if( Lsingularity )
               basis = T{l}{1}(im(j)+1);
            else
               basis = T{l}{1};
            end

            for it=1:nt
              for iz=1:nz
               cosa = cos(im(j)*tarr(it)-in(j)*zarr(iz));
               sina = sin(im(j)*tarr(it)-in(j)*zarr(iz));
               At(:,it,iz) = At(:,it,iz) + basis.*( Ate(l,j)*cosa + Ato(l,j)*sina );
               Az(:,it,iz) = Az(:,it,iz) + basis.*( Aze(l,j)*cosa + Azo(l,j)*sina );
              end
            end
        end
    end


    Acov{1} = At;
    Acov{2} = Az; 
end
