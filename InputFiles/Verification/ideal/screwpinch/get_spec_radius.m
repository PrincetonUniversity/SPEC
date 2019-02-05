function r_out = get_spec_radius(filename, theta, zeta, vol)
% Return the radial position of a KAM surface for a given theta, zeta and
% Nvol
%
% INPUT
% -----
%   filename    hdf5 filename
%   theta:      Poloidal angle
%   zeta:       Toroidal angle
%   vol:        Volume number
%
% OUPUT
% -----
%   r_out:      Radial position of the KAM surface


% Load a bunch of stuff
   
mn     = h5read(filename,'/mn');
im     = h5read(filename,'/im');
in     = h5read(filename,'/in');
Rmn    = h5read(filename,'/Rbc');

r_out = 0;

for k=1:mn
   r_out = r_out + Rmn(k, vol+1) * cos(double(im(k)) * theta - double(in(k)) * zeta);
end

end