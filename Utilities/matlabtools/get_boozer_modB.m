%% get_boozer_modB( BDATA, NTHETA, NPHI  )
% =============================================
% 
% Generates an output data struct for the fourier 
% coefficients Bmn of mod(B) in Boozer coordinates
% 
% INPUT
% -----
%   -bdata    : must be produced by calling read_boozer(filename,root)
%   -Ntheta   : number of meshpoints for theta array
%   -Nphi     : number of meshpoints for phi array
%
% ------------------------------------%
% Written by S.Guinchard (05/18/22)   %
% ------------------------------------%

function modB = get_boozer_modB(b, Ntheta, Nphi)

    xm = double(b.Booz_xForms.Outputs.xm_b);
    xn = double(b.Booz_xForms.Outputs.xn_b);
    BmnB = b.Booz_xForms.Outputs.bmnc_b;
    Nfp  = double(b.Booz_xForms.Inputs.nfp);
    
    theta = linspace(0, 2*pi, Ntheta);
    %phi   = linspace(0, 2*pi/double(Nfp), Nphi); % (for 2D plot of modB boozer)
    phi   = linspace(0, 2*pi, Nphi);             % (for 3D plot of modB boozer)
    
    [modB.Phi, modB.Theta] = meshgrid(phi, theta);
   
    
    modB.modB = zeros(size(modB.Theta));
    for i =1:length(xm)
        m = xm(i);
        n = xn(i);
        angle = (m* modB.Theta - n*modB.Phi);
        modB.modB = modB.modB + BmnB(i)*cos(angle);
            % case of stellarator symmetry (see init.from.spec.py)
            if b.Booz_xForms.Inputs.asym == 1.0
                modB.modB = modB.modB + BmnB(i)*sin(angle);
            end
    end

    
end
   