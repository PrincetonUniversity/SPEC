function B = get_full_field(filename, r, theta, zeta, nr, a)

% Return SPEC magnetic field solution componants as a function of r on the
% entire cylinder. Works in cylindrical geometry.
%
% INPUT
% -----
%   filename:   SPEC .h5 output file
%   r:          r coordinate on which the solution will be interpolated
%   theta:      Theta coordinate 
%   zeta:       Zeta coordinate
%   nr:         Number of points in each volume for interpolation
%   a:          Cylindre radius
%
% OUTPUT
% ------
%   B:          3xlength(r) array containing the B field solution
%               interpolated at r
%
% Written by A.Baillod(2019)

% Load data
fdata = read_spec_field(filename);
Nvol = fdata.Nvol;

epsilon = 0;                   % Minimum radius
B_temp = zeros(3, nr*Nvol);    % Allocate memory
r_temp = zeros(1, nr*Nvol);    % Allocate memory

iimin = 1;
iimax = 0;
for i=1:Nvol
   
   if i==1
       rmin = epsilon;
   else
       rmin = get_spec_radius(filename, 0, 0, i-1);
   end
   
   if i==Nvol
       rmax = a;
   else
       rmax = get_spec_radius(filename, 0, 0, i);
   end
   
   r_vol = linspace(rmin, rmax, nr+1);
   r_vol = r_vol(2:end);
   
   if i==1
        sarr = 2 * ((r_vol - rmin) ./ (rmax - rmin)).^2 - 1;
   else
        sarr = 2 * (r_vol - rmin) ./ (rmax - rmin) - 1;
   end
   
%     Get magnetic field
   temp = get_spec_magfield_cyl(fdata, i, sarr, theta, zeta);
   Bcontrav = zeros(3, length(temp{1}));
   for jj=1:3
      Bcontrav(jj,:) = temp{jj}; 
   end
   
   iimax = iimax + length(r_vol);
   r_temp(iimin:iimax) = r_vol;
   
%    Transform in canonical basis and save in vectors
   B_temp(:, iimin:iimax) = contra2cov_cyl(filename, i, sarr, Bcontrav, nr, ...
                                           theta, zeta, 1);


   iimin = iimin + length(r_vol);
       
end

B = zeros(3,length(r));
for i=1:3
    B(i,:) = interp1(r_temp, B_temp(i,:), r, 'spline');
end





