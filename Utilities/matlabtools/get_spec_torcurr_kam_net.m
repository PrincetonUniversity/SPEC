function Itor = get_spec_torcurr_kam_net(data, ntheta)

%
% GET_SPEC_TORCURR_KAM_NET( DATA, NTHETA )
% ========================================
%
% Calculates the net toroidal surface-current on each interface
%
% INPUT
% -----
%   - data        : obtained via read_spec(filename)
%   - ntheta      : poloidal resolution for the loop integral
%
% OUTPUT
% ------
%   - Itor        : array of size Nvol-1 with values of net-surface-current on each interface (units are mu0*I)
%
% Note: Stellarator symmetry assumed in some of the routines used.
%
% written by J.Loizu (2017)


Mvol     = data.output.Mvol;
Itor     = zeros(1,Mvol-1);

zarr     = 0;
tarr     = linspace(0,2*pi,ntheta);

sarr     = [1 -1];
intB     = [0  0];

for ikam=1:Mvol-1

 lvol     = [ikam  ikam+1];

 for i=1:2

  Bcontrav = get_spec_magfield(data,lvol(i),sarr(i),tarr,zarr);
  gmat     = get_spec_metric(  data,lvol(i),sarr(i),tarr,zarr);

  Bs       = Bcontrav{1};
  Bt       = Bcontrav{2};
  Bz       = Bcontrav{3};

  gst      = gmat{1,2};
  gtt      = gmat{2,2};
  gzt      = gmat{3,2};

  intB(i)  = trapz(tarr, Bs.*gst + Bt.*gtt + Bz.*gzt);

 end

 Itor(ikam) = intB(2)-intB(1);
 
end
