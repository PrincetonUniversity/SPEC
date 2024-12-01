%% get_boozer_coordinates( BDATA, DATA, NPOINTS, NPERIODS )
% =============================================
% 
% Computes the Boozer coordinates from a .boz.h5 file
% .boz.h5 file must be generated using run_boz.py (see corresponding repo)
% 
% INPUT
% -----
%   -bdata : must be produced by calling read_boozer(filename, root)
%   -data  : must be produced using read_spec(filename)
%   -Npoints : number of points to be used for the grid
%   -Nperiods : Number of toroidal periods
%
% OUTPUT
% ------
%   -coordb : matlab struct with 2 fields - phi_boz, theta_boz
% ------------------------------------%
% Written by S.Guinchard (03/29/22)   %
% ------------------------------------%

function coordb = get_boozer_coordinates(b,d, Npoints, Nperiods)

m   = b.Booz_xForms.Outputs.xm_b;
n   = b.Booz_xForms.Outputs.xn_b; 

coord = get_spec_straight_fieldlines(d,Npoints,Nperiods);
iota  = b.Booz_xForms.Inputs.iota;
phi   = coord.phi;
theta = coord.theta_sfl;
numns = b.Booz_xForms.Outputs.numns_b;
theta_b = theta;

for i = 1:Npoints
    for j =1:length(n)
        phi(i) = phi(i) + numns(j) * sin(double(m(j))*theta(i) - double(n(j))*phi(i));
        theta_b(i) = theta_b(i) + iota * numns(j) * sin(double(m(j))*theta(i) - double(n(j))*phi(i));
    end
end

coordb.phib   = wrapTo2Pi(phi);
coordb.thetab = wrapTo2Pi(theta_b);

end