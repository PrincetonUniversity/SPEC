function B = get_full_field(data, r, theta, zeta, nr)

%
% GET_FULL_FIELD( DATA, R, THETA, ZETA, NR )
% ==========================================
%
% Return SPEC magnetic field solution componants as a function of r, theta,
% zeta. The coordinate r is constructed from the radial position of the
% volume interface for a given pair (theta, zeta) and from the coordinate
% s.
%
% INPUT
% -----
%   data:   	data obtained from read_spec(filename)
%   r:          Radial coordinate (array of doubles)
%   theta:      Theta coordinate (double)
%   zeta:       Zeta coordinate (double)
%   nr:         Number of points in each volume for radial interpolation
%
% OUTPUT
% ------
%   B:          3xlength(r)xlength(theta)xlength array containing the B field solution
%               interpolated at r
%
% Written by A.Baillod(2019)


% Load data
G = data.input.physics.Igeometry;
Nvol = data.output.Mvol;

epsilon = 1E-16;

nt = length(theta);
nz = length(zeta);

B_temp = zeros(3, nr*Nvol, nt, nz);    % Allocate memory
r_temp = zeros(1, nr*Nvol, nt, nz);    % Allocate memory

iimin = 1;
iimax = 0;

r0 = get_spec_radius(data, theta, zeta, 0);
B = zeros(3,length(r),nt,nz);
for i=1:Nvol
   
   % if first volume, don't take rmin=0
   if i==1
       rmin = epsilon;
   else
       [ri, zi] = get_spec_radius(data, theta, zeta, i-1);
       rmin = sqrt((ri-r0)^2+zi^2);
   end
   
   [ri, zi] = get_spec_radius(data, theta, zeta, i);
   rmax = sqrt((ri-r0)^2+zi^2);
   
   r_vol = linspace(rmin, rmax, nr+1);  % Minor radius
   r_vol = r_vol(2:end);
   
   if i==1
       if G == 1
        sarr = 2.0 * (r_vol - rmin) ./ (rmax - rmin) - 1; 
       else
        sarr = 2.0 * ((r_vol - rmin) ./ (rmax - rmin)).^2 - 1;
       end
   else
       sarr = 2.0 * (r_vol - rmin) ./ (rmax - rmin) - 1;
   end
   Bcontrav = get_spec_magfield(data, i, sarr, theta, zeta);
   
   % Generate radial coordinate array
   iimax = iimax + length(r_vol);
   r_temp(iimin:iimax) = r_vol;
   
   % And convert it to covariant basis (normalized here)
   B_cov = contra2cov(data, i, Bcontrav, sarr, theta, zeta, 1);
   
   B_temp(1, iimin:iimax, :, :) = B_cov{1};
   B_temp(2, iimin:iimax, :, :) = B_cov{2};
   B_temp(3, iimin:iimax, :, :) = B_cov{3};
                                    
   iimin = iimin + length(r_vol); 
                    
   
   for comp = 1:3
       ind = find(r<=rmax);
       r_int = r(ind);
       ind = find(r_int>=rmin);
       r_int = r_int(ind);
       
       for j = 1:nt
           for k = 1:nz
            B(comp, ind, j, k) = interp1(r_vol, B_cov{comp}(:,j,k), r_int, 'spline');
           end
       end
   end
end






