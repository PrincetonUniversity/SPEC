function [sarr, R] = getR_derivatives(filename, vol, ns, theta, phi)
%
% Returns the derivatives of R corresponding to SPEC coordinate 
% system
%
% INPUT
% -----
%	filename: 	SPEC hdf5 output filename
%	vol:		Volume number
%	ns:		s-coordinate resolution
%	theta:		Theta angle
%	phi:		Phi angle
%
% OUTPUT
% ------
%	sarr:		s-coordinate array
%	R:		4xlength(sarr) array containing R, dR / ds, 
%			dR / dtheta and dR / dphi
%
% Written by A.Baillod (2019)


% Load geometry
mn     = h5read(filename,'/mn');
im     = h5read(filename,'/im');
in     = h5read(filename,'/in');
Rmn = h5read(filename,'/Rbc');

% Allocate data for R and its derivative in s, theta and phi (4), for each
% and for ns points 
Rarr = zeros(4, ns); 

% s array declaration
sarr = linspace(-1, 1, ns);
sbar = (sarr + 1) / 2.0;
                        
                        
for imn=1:mn
    cosa = cos(double(im(imn)*theta - in(imn)*phi));
    sina = sin(double(im(imn)*theta - in(imn)*phi));
    
    if vol==1
        if im(imn)==0
            factor  = sqrt(sbar);
            dfactor  = 0.25./sqrt(sbar);
        else
            factor = sbar.*(double(im(imn)) / 2); 
            dfactor  = (double(im(imn))/4)*sbar.^(double(im(imn))/2-1);
        end
    else
        factor  = sbar;
        dfactor  = 1.0 / 2.0 *ones(1,length(sbar));
    end


    Rarr(1,:) = Rarr(1,:)   + (Rmn(imn, vol) + (Rmn(imn, vol+1) ...
                            - Rmn(imn, vol))*factor) * cosa;
    Rarr(2,:) = Rarr(2,:) + (Rmn(imn, vol+1) ...
                          - Rmn(imn, vol)) * dfactor * cosa;
    Rarr(3,:) = Rarr(3,:) - (Rmn(imn, vol) + (Rmn(imn, vol+1) ...
                          - Rmn(imn, vol))*factor) * double(im(imn)) *sina;
    Rarr(4,:) = Rarr(4,:) + (Rmn(imn, vol) + (Rmn(imn, vol+1) ...
                          - Rmn(imn, vol))*factor) * double(in(imn)) *sina;
end

R = Rarr;
