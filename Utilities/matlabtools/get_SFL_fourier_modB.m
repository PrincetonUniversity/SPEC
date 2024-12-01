%% get_SFL_fourier_modB( DATA, NTHETA, NPHI,  )
% =============================================
% 
% Generates an output data struct for the fourier 
% coefficients Bmn of mod(B)in SFL coordinates
% 
% INPUT
% -----
%   -data     : must be produced by calling read_spec(filename)
%   -Ntheta   : number of meshpoints for theta array
%   -Nphi     : number of meshpoints for phi array
%
% ------------------------------------%
% Written by S.Guinchard (05/23/22)   %
% ------------------------------------%
function data_Bmn = get_SFL_fourier_modB(d, Ntheta, Nphi)


    Nfp  = d.input.physics.Nfp;
    Ntor = d.input.physics.Ntor;
    Mpol = d.input.physics.Mpol;
        K = 2*pi^2/double(Nfp);

    lambdamn = squeeze(d.output.lambdamn);
    lambdamn = lambdamn(1:end,2);
    mprime   = d.output.ims;
    nprime   = d.output.ins/double(Nfp);

    data_Bmn.Theta = linspace(0, 2*pi, Ntheta);
    %data_Bmn.Phi   = linspace(0, 2*pi/double(Nfp), Nphi); %change accordingly 
    data_Bmn.Phi   = linspace(0, 2*pi, Nphi);


    data_Bmn.modB = zeros(Ntheta, Nphi);
    modB_SPEC = squeeze(get_spec_modB(d,1,1,data_Bmn.Theta,data_Bmn.Phi));

    sumlmn = 0.0;

    for ind =2:length(lambdamn)
        for line = 1:Ntheta
            for column = 1:Nphi
                sumlmn = sumlmn + mprime(ind)*lambdamn(ind)*cos(double(mprime(ind))*data_Bmn.Theta(line)-double(nprime(ind))*data_Bmn.Phi(column));
            end
        end
    end

    data_Bmn.jacobian = abs(double(1+sumlmn));

    for m=0:Mpol
        for n=-Ntor:Ntor
            if m == 0 && n<0
               continue 
            end
            for line = 1:Ntheta
                for column = 1:Nphi
                    % evaluate modB_SPEC*cos(m theta - n phi) on the grid (theta, phi)
                    data_Bmn.modBcos(line,column) = modB_SPEC(line,column)* data_Bmn.jacobian ...
                                * cos(double(m)*data_Bmn.Theta(line) - double(Nfp)*double(n)*data_Bmn.Phi(column));
                 end
            end
            
            data_Bmn.Bmn(m+1,n+Ntor+1) = 1/K*trapz(data_Bmn.Phi,trapz(data_Bmn.Theta,data_Bmn.modBcos,1));
        end
    end
% introduce missing 1/2 coefficient for mode (0,0)
    data_Bmn.Bmn(1,Ntor+1) = 1/2*data_Bmn.Bmn(1,Ntor+1);
    
    for i=1:Ntheta
        for m=0:Mpol
            for n=-Ntor:Ntor
                if m == 0 && n<0
                    continue 
                end
                data_Bmn.modB(i,:) = data_Bmn.modB(i,:) + data_Bmn.Bmn(m+1,n+Ntor+1)* ...
                                     cos(double(m)*data_Bmn.Theta(i) - double(Nfp)*double(n)*data_Bmn.Phi);
            end
        end
    end
    data_Bmn.modB = data_Bmn.modB/2;

%end of function
end