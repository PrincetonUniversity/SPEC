function [r_out, z_out] = get_spec_radius(data, theta, zeta, vol)

%
% GET_SPEC_RADIUS( DATA, THETA, ZETA, VOL )
% =========================================
%
% Return the radial position of a KAM surface for a given theta, zeta and
% Nvol
%
% INPUT
% -----
%   data        Obtained via read_spec(filename)
%   theta:      Poloidal angle
%   zeta:       Toroidal angle
%   vol:        Volume number
%
% OUPUT
% -----
%   r_out:      Radial position of the KAM surface


% Load a bunch of stuff
   
mn     = data.output.mn;
im     = data.output.im;
in     = data.output.in;
Rmn    = data.output.Rbc;
Zmn    = data.output.Zbs;
G      = data.input.physics.Igeometry;


switch G
    case 1
        r_out = 0;

        for k=1:mn
           r_out = r_out + Rmn(k, vol+1) * cos(double(im(k)) * theta - double(in(k)) * zeta);
        end
        z_out = 0;
    case 2
        r_out = 0;

        for k=1:mn
           r_out = r_out + Rmn(k, vol+1) * cos(double(im(k)) * theta - double(in(k)) * zeta);
        end
        z_out = 0;
    case 3
        r_out = 0;
        z_out = 0;
        
        for k=1:mn
           r_out = r_out + Rmn(k, vol+1) * cos(double(im(k)) * theta - double(in(k)) * zeta);
           z_out = z_out + Zmn(k, vol+1) * sin(double(im(k)) * theta - double(in(k)) * zeta);
        end

end
