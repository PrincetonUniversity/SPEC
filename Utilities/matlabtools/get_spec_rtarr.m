function rtdata = get_spec_rtarr(data,lvol,sarr,tarr,zarr0)
 
 
% Transforms (s,theta) array into (R,theta) array in volume number lvol in slab or cylindrical geometry
%
% INPUT
%   -data    : must be produced by calling e.g. read_spec_grid(filename)
%   -lvol    : volume number
%   -sarr    : is the array of values for the s-coordinate
%   -tarr    : is the array of values for the theta-coordinate
%   -zarr    : is the array of values for the zeta-coordinate
%
% OUTPUT
%   -rtdata  : array with (R,theta,dRds) data array with size 3*ns*nt where ns=length(sarr),nt=length(tarr)
%
% Note: Stellarator symmetry is assumed
%
% written by J.Loizu (2018)
% updated by J.Loizu (2020)


Rac     = data.Rbc(:,lvol);   % inner volume boundary harmonics
Rbc     = data.Rbc(:,lvol+1); % outer volume boundary harmonics


sarr    = transpose(sarr);
ns      = length(sarr);
nt      = length(tarr);
sbar    = (sarr+1)/2;

mn      = data.mn;
im      = double(data.im);
in      = double(data.in);

Rarr    = zeros(ns,nt); % allocate data for R-array
Tarr    = zeros(ns,nt); % allocate data for theta-array
dRarr   = zeros(ns,nt); % allocate data for R-array derivative (in s)

fac     = cell(mn,2);   % allocate data for regularization factors 


% Construct regularization factors

switch data.Igeometry
 case 1
  for j=1:mn
   fac{j}{1} = sbar;
   fac{j}{2} = 0.5;
  end
 case 2
  for j=1:mn
    if lvol==1
        if im(j)==0
            fac{j}{1}  = sqrt(sbar);
            fac{j}{2}  = 0.25./sqrt(sbar);
        else
            fac{j}{1} = sbar.*(im(j) / 2);
            fac{j}{2}  = (im(j)/4)*sbar.^(im(j)/2-1);
        end
    else
        fac{j}{1}  = sbar;
        fac{j}{2}  = 0.5;
    end
  end
end


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