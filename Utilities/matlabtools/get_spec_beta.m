function [beta_ax, beta_av] = get_spec_beta(data, vols)

%
% GET_SPEC_BETA( DATA, VOLS )
% ===========================
%
% Calculates beta of the equilibrium, both the average beta=2*<p>/B(0)^2 
% and the axis beta=2*p(0)/B(0)^2, and returns the latter
%
% INPUT
% -----
%   -data    : Produced via read_spec('filename.sp.h5');
%   -vols    : volumes on which beta is averaged
%
% OUTPUT
% ------
%  -beta_ax  : value of beta on axis
%  -beta_av  : value of beta on average
%
% written by J.Loizu (2016) 
% modified by J.Loizu (05.2017)
% modified by A.Baillod (06.2019)
%
% A.Baillod: update for new hdf5 format (04.2020)
%

% Read some data
pscale = data.input.physics.pscale;
press  = pscale * data.input.physics.pressure;

% Define integration coordinates
tarr    = linspace(0,2*pi,16);
zarr    = linspace(0,2*pi,16);

% Allocate for average beta values
av_beta = zeros(1,length(vols));
volume =0;

% In each volume, get the volume, the jacobian and the field modulus, and
% integrate J/B^2. Store each beta averaged on a volume in beta_av
for lvol=1:length(vols)
% disp(['Volume ', num2str(lvol), ' out of ', num2str(length(vols)), '.'])
 volume_number = vols(lvol);
 
 if volume_number==1
    sarr    = linspace(-0.99,  1, 5);
 else
    sarr    = linspace(-1   ,  1, 5);
 end
     
 % Compute total volume
 volume = volume + get_spec_volume(data,volume_number,64,64,64);
 
 % And integral 2*p/B^2
 modB     = get_spec_modB(data,volume_number,sarr,tarr,zarr);
 jacobian = get_spec_jacobian(data, volume_number, sarr, tarr, zarr);
 arg = jacobian ./ (modB.^2);
 av_beta(lvol) = 2*press(volume_number)*trapz(zarr, trapz(tarr, trapz(sarr, arg, 1), 2), 3);
 
end
beta_av = sum(av_beta) / volume;

ind = find(vols==1);
if length(ind)==0
    volume_number = 1;
    volume = get_spec_volume(data,volume_number,64,64,64);
    modB     = get_spec_modB(data,volume_number,sarr,tarr,zarr);
    jacobian = get_spec_jacobian(data, volume_number, sarr, tarr, zarr);
    arg = jacobian ./ (modB.^2);
    beta_ax = av_beta(1);
else
    beta_ax = av_beta(ind);
end
