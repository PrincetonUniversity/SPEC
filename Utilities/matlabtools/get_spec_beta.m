function beta = get_spec_beta(fname, vol_ind)
 
% Calculates beta of the equilibrium, both the average beta=2*<p>/B(0)^2 
% and the axis beta=2*p(0)/B(0)^2, and returns the latter
%
% INPUT
%  -fname   : path to the hdf5 output file (e.g. 'testcase.sp.h5')
%
% OUTPUT
%  -beta		: value of beta on axis
%
% written by J.Loizu (2016) 
% modified by J.Loizu (05.2017)


gdata  = read_spec_grid(fname);

fdata  = read_spec_field(fname);

Nvol   = h5read(fname,'/Nvol');

pscale = h5read(fname,'/pscale');

volume = zeros(Nvol,1);

press  = pscale*h5read(fname,'/pressure');

for lvol=1:length(vol_ind)
    
 ivol = vol_ind(lvol);

 volume(lvol) = get_spec_volume(gdata,ivol,64,64,64);

end

avpress = sum(press.*volume)/sum(volume);

zarr    = linspace(0,2*pi,64);
modB    = get_spec_modB(fdata,1,-0.99,0,zarr);
B0      = mean(modB);

beta_ax = 2*press(1)/(B0^2)

beta_av = 2*avpress/(B0^2)

beta 		= beta_ax;

