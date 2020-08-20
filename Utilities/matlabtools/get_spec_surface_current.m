function [tflux, IPDt] = get_spec_surface_current(data, ns, nt, zeta)

% Returns the sheet current flowing through each interface, normalized by 
% mu_0. This routine computes the actual integral of the poloidal field to 
% compute the current via Ampere's law. This requires to specify zeta; 
% however, since the surface current is a flux function, zeta should not 
% have any influence on the result.
%
% INPUT
% -----
%   data:   	data obtained via read_spec(filename)
%   ns:         Radial resolution
%   nt:         Poloidal resolution
%   zeta:       Toroidal angle - see introduction remark
%
% OUTPUT
% ------
%   tflux:      The toroidal flux enclosed by each interface (1xNvol)
%   IPDt:       The toroidal surface current in each interface (1xNvol-1),
%               normalized by mu_0
%
% Written by A.Baillod (2019)


% Constant definition
mu0 = 4*pi*1E-7;
epsilon = 1E-5;

% Data loading
fdata = fdata_from_data(data);      % Read data
Nvol = fdata.Mvol;                      % Total number of volumes
sarr = linspace(-1, 1, ns);

% Allocate memory
Bcov = cell(1, Nvol);

theta = linspace(0, 2*pi, nt);

IPDt = zeros(1, Nvol-1);

% Get magnetic field
for ivol=1:Nvol
    if ivol==1
        sarr(1)=-1+epsilon;
    else
        sarr(1)=-1;
    end
    
    temp = get_spec_magfield(fdata, ivol, sarr, theta, zeta);  
    
    Bcov{ivol} = contra2cov(fdata, ivol, temp, sarr, theta, zeta, 0);                     
end

for ivol=1:Nvol-1
    dBtheta = -Bcov{ivol}{2}(end,:,1) + Bcov{ivol+1}{2}(1,:,1);    
    IPDt(ivol) = trapz(theta, dBtheta);
end

tflux = fdata.tflux;

end
