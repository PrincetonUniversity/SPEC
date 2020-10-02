function [psi_coord, I_vol] = get_spec_volume_current(data, cumul)

%
% GET_SPEC_VOLUME_CURRENT( DATA, CUMUL )
% ======================================
%
% Returns the toroidal current flowing in each volume, normalized by mu_0
%
% INPUT
% -----
%   data:       data obtained via read_spec(filename)
%
% OUTPUT
% ------
%   psi_coord:  The toroidal flux enclosed by each interface
%   I_vol:      The toroidal volume current flowing in each volume
%               (cumulative)
%
% Written by A.Baillod (2019)


% Data loading
Nvol = data.input.physics.Nvol;                      % Total number of volumes

% Data processing

% First, get the current in each volume
psi_coord = zeros(1, Nvol);             % Allocate memory
I_vol = zeros(1, Nvol);

mu = data.output.mu;
tflux = data.output.tflux;
sumI = 0;
phiedge = data.input.physics.phiedge;
    
for ivol=1:Nvol

    if ivol==1    
        I_vol(ivol) = mu(ivol) * tflux(ivol) * phiedge;
    else
        % Add previous current volumes (sumI) since we use a cumulative 
        % representation
        if cumul
            I_vol(ivol) = mu(ivol)  * (tflux(ivol) - tflux(ivol-1)) * phiedge + sumI;
        else
            I_vol(ivol) = mu(ivol)  * (tflux(ivol) - tflux(ivol-1)) * phiedge;
        end
    end
    
    psi_coord(ivol) = tflux(ivol);    
    
    sumI = I_vol(ivol);
end

end 
