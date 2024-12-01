function W = get_spec_energy_slab(data,lv,ns,nt,nz)
 
%
% GET_SPEC_ENERGY_SLAB( DATA, LV, NS, NT, NZ )
% ============================================
%
% Calculates Plasma Magnetic Energy in volume lv in slab geometry 
%
% INPUT
% -----
%   -data      : must be produced by calling read_spec(filename)
%   -lv        : volume in which to calculate the energy (total energy for lvol=0)
%   -ns        : is the s-coordinate resolution 
%   -nt        : is the theta-coordinate resolution
%   -nz        : is the zeta-coordinate resolution
%
% OUTPUT
% ------
%   -W         : total energy 
%
% Note: Stellarator symmetry is assumed
%
% written by J.Loizu (2018)

    % Test input
    Igeometry = data.input.physics.Igeometry;
    if Igeometry~=1
        error('This routine is only for slab geometries')
    end

    % Read data
    Nvol = data.input.physics.Nvol;

    W    = 0;

    sarr = linspace(-1,1,ns);
    tarr = linspace(0,2*pi,nt);
    zarr = linspace(0,2*pi,nz);

    if(lv==0)
        lvol_start=1;
        lvol_end  =Nvol;
    else
        lvol_start=lv;
        lvol_end  =lv;
    end
    
    
    for lvol=lvol_start:lvol_end

        % Evaluate mod B and jacobian
        modB = get_spec_modB(data,lvol,sarr,tarr,zarr);
        jac  = get_spec_jacobian(data,lvol,sarr,tarr,zarr);

        % Integrate
        F    = jac.*modB.^2;

        dW   = trapz(sarr,trapz(tarr,trapz(zarr,F,3),2));
        W    = W + 0.5*dW;
    end


