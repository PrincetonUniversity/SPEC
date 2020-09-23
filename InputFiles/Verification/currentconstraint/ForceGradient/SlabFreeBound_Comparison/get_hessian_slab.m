function hessian = get_hessian_slab(Nvol, tflux, pflux, mu, R, Isurf, Icoil, debug)

% fname = 'Slab_FreeBound_Nvol1.sp.h5';
% data = read_spec(fname);
% 
% R = data.output.Rbc;
% Nvol = data.input.physics.Nvol;
% tflux = data.output.tflux; tflux(2:end) = tflux(2:end) - tflux(1:end-1);
% pflux = data.output.pflux;
% mu = data.output.mu;
% Isurf = data.input.physics.Isurf;
% Ivolume = data.input.physics.Ivolume; Ivolume(2:end) = Ivolume(2:end) - Ivolume(1:end-1);
% Icoil = data.input.physics.curpol;


switch Nvol
    case 1
        cosmu = cos(mu(1)*R(2));
        sinmu = sin(mu(1)*R(2));
        
        M(1,1) =  cosmu - 1;
        M(2,1) =  sinmu    ;
        M(1,2) =  sinmu    ;
        M(2,2) = -cosmu + 1;
        
        RHS(1) =  mu(1) * tflux(1) / (2*pi);
        RHS(2) =  mu(1) * pflux(1) / (2*pi);
        
        x      = M\RHS';
        A      = x(1);
        B      = x(2);
        
        BPthe  = A*cosmu + B*sinmu;
        BPphi  = B*cosmu - A*sinmu;
        BP     = sqrt(A^2 + B^2);
        
        BVthe  = Isurf(1) / (2*pi) + BPthe;
        BVphi  = Icoil / (2*pi);
        BV     = sqrt(BVthe^2 + BVphi^2);
        
        force  = .5 * (BV^2 - BP^2);        %OK, VERIFIED, matched force in SPEC
        
        dM(1,1) = -sinmu * mu(1);
        dM(1,2) =  cosmu * mu(1);
        dM(2,1) =  cosmu * mu(1);
        dM(2,2) =  sinmu * mu(1);
        
        dx = - M\(dM * x);
        dA = dx(1);
        dB = dx(2);
        
        dBBPdA = 2*A;
        dBBPdB = 2*B;
        
        dBVphi = 0;
        dBVthe = dA * cosmu - A * mu(1) * sinmu + dB * sinmu + B * mu(1) * cosmu;
        
        dBBP = dBBPdA * dA + dBBPdB * dB;                   % OK, matches dFFdRZ(1,1,1,1,1) *2;
        dBBV = 2 * BVthe * dBVthe + 2 * BVphi * dBVphi;     % OK
        
        hessian = .5 * (dBBV - dBBP);
        

    otherwise
        
        disp('Unsupported number of volumes')
end


    
    
    
    
end


