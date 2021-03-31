function rtdata = get_spec_rtarr(data,lvol,sarr,tarr,zarr0)
 
%
% GET_SPEC_RTARR( DATA, LVOL, SARR, TARR, ZARR0 )
% ===============================================
%
% Transforms (s,theta) array into (R,theta) array in volume number lvol in slab or cylindrical geometry
%
% INPUT
% -----
%   -data    : must be produced by calling e.g. read_spec(filename)
%   -lvol    : volume number
%   -sarr    : is the array of values for the s-coordinate
%   -tarr    : is the array of values for the theta-coordinate
%   -zarr    : is the array of values for the zeta-coordinate
%
% OUTPUT
% ------
%   -rtdata  : array with (R,theta,dRds) data array with size 3*ns*nt where ns=length(sarr),nt=length(tarr)
%
% Note: Stellarator symmetry is assumed
%
% written by J.Loizu (2018)
% updated by J.Loizu (2020)


Rac     = data.output.Rbc(:,lvol);   % inner volume boundary harmonics
Rbc     = data.output.Rbc(:,lvol+1); % outer volume boundary harmonics

if(size(sarr,1)==1)
sarr    = transpose(sarr);
end

ns      = length(sarr);
nt      = length(tarr);
sbar    = (sarr+1)/2;

mn      = data.output.mn;
im      = double(data.output.im);
in      = double(data.output.in);

Rarr    = zeros(ns,nt); % allocate data for R-array
Tarr    = zeros(ns,nt); % allocate data for theta-array
dRarr   = zeros(ns,nt); % allocate data for R-array derivative (in s)



% Construct regularization factors

fac = get_spec_regularization_factor(data, lvol, sarr, 'G');

% Construct (R,theta) coordinates array

for j=1:mn
  for it=1:nt
     cosa = cos(im(j)*tarr(it)-in(j)*zarr0);
     sina = sin(im(j)*tarr(it)-in(j)*zarr0);
     Rarr(:,it)  = Rarr(:,it) + (Rac(j) + fac{j}{1}.*(Rbc(j)-Rac(j)) )*cosa;
     dRarr(:,it) = dRarr(:,it) + fac{j}{2}*(Rbc(j)-Rac(j))*cosa;
     Tarr(:,it)  = tarr(it);
  end
end

rtdata{1} = Rarr;
rtdata{2} = Tarr;
rtdata{3} = dRarr;
