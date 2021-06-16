function vcov = contra2cov(data, vol, vcontrav, s, theta, phi, norm)

%
% CONTRA2COV( DATA, VOL, VCONTRAV, S, THETA, PHI, NORM )
% ======================================================
%
% Transform a contravariant vector to a covariant one in any geometry
%
% INPUT
% -----
%   data:       via read_spec(filename)
%   vol:        volume
%   vcontrav:   structure containing contravariant components as a function of r, size = 3xlength(r)
%   s:          radial coordinate
%   theta:      theta angle
%   phi:        phi angle
%   norm:       use the unitary (canonical) basis (=1) or the general 
%               basis(=0)
%
% OUTPUT
% ------
%   vcov:       covariant vector as a function v

% Read geometry
G = data.input.physics.Igeometry;


% Get R derivatives
[s, R] = get_spec_R_derivatives(data, vol, s, theta, phi, 'R');
[s, Z] = get_spec_R_derivatives(data, vol, s, theta, phi, 'Z');


nt = length(theta);
np = length(phi);
ns = length(s);

% transform in covariant basis    
g = get_spec_metric(data, vol, s, theta, phi);
        
switch G
    case 1
        norm1 = sqrt(1 + R{3}.^2 + R{4}.^2) ./ R{2};
        norm2 = ones(size(R{1}));
        norm3 = ones(size(R{1}));
        
    case 2                     
        norm1 = sqrt(R{3}.^2 + R{1}.^2 + (R{1}.*R{4}).^2) ./ (R{2}.*R{1});
        norm2 = 1 ./ R{1};
        norm3 = ones(size(R{1}));
        
    case 3        
        norm1 = sqrt((R{3}.*Z{4}).^2 + (R{4}.*Z{3}).^2 + (R{1}.*Z{3}).^2  ...
               -2*R{3}.*R{4}.*Z{3}.*Z{4} + (R{1}.*R{3}).^2);
        norm2 = sqrt((R{4}.*Z{2}).^2 + (R{1}.*Z{2}).^2 + (R{2}.*Z{4}).^2 + (R{1}.*R{2}).^2 ...
               -2*R{2}.*R{4}.*Z{2}.*Z{4});
        norm3 = sqrt((R{2}.*Z{3}).^2 + (R{3}.*Z{2}).^2 -2*R{2}.*R{3}.*Z{2}.*Z{3});
        
end

% Compute covariant components
temp = cell(1,3);
vcov = cell(1,3);

for ii=1:3
   temp{ii} = zeros(ns, nt, np); 
   vcov{ii} = zeros(ns, nt, np);
end



for is = 1:ns
    for it = 1:nt
        for ip = 1:np
            temp{1} = vcontrav{1}.*g{1}{1} + vcontrav{2}.*g{1}{2} + vcontrav{3}.*g{1}{3};
            temp{2} = vcontrav{1}.*g{2}{1} + vcontrav{2}.*g{2}{2} + vcontrav{3}.*g{2}{3};
            temp{3} = vcontrav{1}.*g{3}{1} + vcontrav{2}.*g{3}{2} + vcontrav{3}.*g{3}{3};
        end
    end
end

% Normalize to canonical basis if necessary
if norm
   vcov{1} = temp{1} .* norm1; 
   vcov{2} = temp{2} .* norm2; 
   vcov{3} = temp{3} .* norm3; 
else
   vcov = temp;
end


end




