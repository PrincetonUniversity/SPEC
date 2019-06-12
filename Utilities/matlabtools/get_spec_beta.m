function [beta_ax, beta_av] = get_spec_beta(fname, vol_ind)
 
% Calculates beta of the equilibrium, both the average beta=2*<p>/B(0)^2 
% and the axis beta=2*p(0)/B(0)^2
%
% INPUT
%  -fname   : path to the hdf5 output file (e.g. 'testcase.sp.h5')
%  -volumes : Volumes in which the average should be done
%
% OUTPUT
%  -beta_ax : value of beta on axis
%  -beta_av : beta average
%
% written by J.Loizu (2016) 
% modified by J.Loizu (05.2017)
% modified by A.Baillod (06.2019)


gdata  = read_spec_grid(fname);
fdata  = read_spec_field(fname);

Nvol   = h5read(fname,'/Nvol');
pscale = h5read(fname,'/pscale');
press  = pscale*h5read(fname,'/pressure');

volume = zeros(Nvol,1);

for lvol=1:length(vol_ind)
    
 ivol = vol_ind(lvol);

 volume(lvol) = get_spec_volume(gdata,ivol,64,64,64);

end

avpress = sum(press(vol_ind).*volume)/sum(volume);

zarr    = linspace(0,2*pi,64);
modB    = get_spec_modB(fdata,1,-0.99,0,zarr);
B0      = mean(modB);

beta_ax = 2*press(1)/(B0^2);

beta_av = 2*avpress/(B0^2);

