%% GET_SPEC_STRAIGHT_FIELDLINES( d, NPOINTS, Nperiods )
% =======================================================
%
% Gives as output the straightfieldlines coordinates from SPEC out file
% 
% INPUT
% -----
%   -data     : must be produced by calling read_boozer(filename, root)
%   -Npoints  : number of points for toroidal resolution
%   -Nperiods : number of toroidal periods
%
% ------------------------------------%
% Written by S.Guinchard (05/17/22)   %
% ------------------------------------%
function SFL_coord = get_spec_straight_fieldlines(d,Npoints,Nperiods)

      m = double(d.output.ims);
      n = double(d.output.ins);
      lambda_mn = d.output.lambdamn(1:end,1,2);
      
      coord = plot_spec_fieldlines(d,Npoints,Nperiods,0);
      theta = coord.theta;
      phi   = coord.phi;
      
      theta_sfl = theta;

        for i = 1:length(phi)
            for j = 1:length(lambda_mn)
                
                theta_sfl(i) = theta_sfl(i) + lambda_mn(j,1,1)*sin(m(j)*theta(i) - n(j)*phi(i));

            end
        end
   
        SFL_coord.theta_sfl = wrapTo2Pi(theta_sfl);
        SFL_coord.phi = wrapTo2Pi(phi);
end