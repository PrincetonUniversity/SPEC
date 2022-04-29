function rmn = get_spec_fourier( data, surface_index, Mpol, Ntor )
%
% GET_SPEC_FOURIER( DATA, SURFACE_INDEX, MPOL, NTOR )
% ===================================================
%
% Extract the fourier harmonics from a SPEC Poincaré data
%
% INPUTS
% ------
%  - DATA: obtained via read_spec()
%  - SURFACE_INDEX: index of surface you want to target
%  - MPOL: Poloidal resolution
%  - NTOR: Toroidal resolution
%
% OUTPUT
% ------
%  - RMN: fluxSurface instance containing the fourier harmonics of the flux
%         surface 
%
% 

       Nfp = double(data.input.physics.Nfp);
       stellsym = 1;
       rmn = fluxSurface( Nfp, Mpol, Ntor, stellsym );

       
       R = data.poincare.R( surface_index, :, : );
       Z = data.poincare.Z( surface_index, :, : );
       
       s = size(R);
       
       R = reshape(R, s(2), s(3));
       Z = reshape(Z, s(2), s(3));
       
       s = size(R);       
       nz = s(1); 
       
       phi  = (0:nz-1.0)*(2.0*pi/nz)/Nfp;
       
       rmn_axis = data.output.Rbc(1:data.input.physics.Ntor+1,1);
       zmn_axis = data.output.Zbs(1:data.input.physics.Ntor+1,1);
    
       r_axis = zeros(1, nz);
       z_axis = zeros(1, nz);
       for n=0:double(data.input.physics.Ntor)
       
           r_axis = r_axis + rmn_axis(n+1) * cos(n*Nfp*phi);
           z_axis = z_axis - zmn_axis(n+1) * sin(n*Nfp*phi);
           
       end
       
       
       rmn = rmn.initialize_from_poincare( R, Z, phi, r_axis, z_axis );

end