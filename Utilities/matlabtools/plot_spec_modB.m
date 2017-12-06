function rzbdata = plot_spec_modB(fname,lvol,sarr,tarr,zarr,newfig)

% Produces plot of |B| in (R,Z,zarr) cross-section(s)
%
% INPUT
%   -fname   : path to the hdf5 output file (e.g. 'testcase.sp.h5')
%   -lvol    : volume number
%   -sarr    : is the array of values for the s-coordinate ('d' for default)
%   -tarr    : is the array of values for the theta-coordinate ('d' for default)
%   -zarr    : is the array of values for the zeta-coordinate ('d' for default)
%   -newfig  : opens(=1) or not(=0) a new figure
%
% OUTPUT
%   -rzbdata : cell structure with 3 arrays: R-data, Z-data, |B|-data
%
% written by J.Loizu (2016)

if(sarr=='d')
sarr=linspace(-1,1,64);
end

if(tarr=='d')
tarr=linspace(0,2*pi,64);
end

if(zarr=='d')
zarr=0;
end

rzbdata = cell(3);


% Read vector potential

fdata  = read_spec_field(fname);

% Compute |B|

modB   = get_spec_modB(fdata,lvol,sarr,tarr,zarr);

% Compute function (R,Z)(s,theta,zeta)

rzdata = get_spec_rzarr(fdata,lvol,sarr,tarr,zarr);

R = rzdata{1};   
Z = rzdata{2};

% Plot

for iz=1:length(zarr)
 if(newfig==1)
 figure
 end
 
 hold on
 
 pcolor(R(:,:,iz),Z(:,:,iz),modB(:,:,iz)); shading interp; colorbar

 axis equal
 title('|B|');
 xlabel('R');
 ylabel('Z');
end

% Output data

rzbdata{1} = R;
rzbdata{2} = Z;
rzbdata{3} = modB;

