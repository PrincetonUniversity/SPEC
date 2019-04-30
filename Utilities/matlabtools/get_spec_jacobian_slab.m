function jacobian = get_spec_jacobian_slab(data,lvol,sarr,tarr,zarr)
 
 
% Calculates Jacobian in volume number lvol for slab geometry
%
% INPUT
%   -data      : must be produced by calling e.g. read_spec_grid(filename)
%   -lvol      : volume number
%   -sarr      : is the array of values for the s-coordinate
%   -tarr      : is the array of values for the theta-coordinate
%   -zarr      : is the array of values for the zeta-coordinate
%
% OUTPUT
%   -jacobian  : Jacobian of the coordinates with size ns*nt*nz where ns=length(sarr),nt=length(zarr),nt=length(zarr)
%
% Note: Stellarator symmetry is assumed
%
% written by J.Loizu (2017)


Rac     = data.Rbc(:,lvol);   % inner volume boundary harmonics
Rbc     = data.Rbc(:,lvol+1); % outer volume boundary harmonics

rpol    = data.rpol; % poloidal extent of the slab is 2*pi*rpol
rtor    = data.rtor; % toroidal extent of the slab is 2*pi*rtor

sarr    = transpose(sarr);
ns      = length(sarr);
nt      = length(tarr);
nz      = length(zarr);
sbar    = (sarr+1)/2;

mn      = data.mn;
im      = double(data.im);
in      = double(data.in);

Rarr    = cell(3); % allocate data for R-array and its derivatives in (s,theta)

for k=1:3
Rarr{k} = zeros(ns,nt,nz); 
end

fac     = cell(mn,2);      % allocate data for regularization factors and their derivatives


% Construct radial factors and their derivatives

for j=1:mn
 fac{j}{1}  = sbar;
 fac{j}{2}  = 0.5*ones(ns,1);
end


% Construct R array and its derivatives

for j=1:mn
  for it=1:nt
    for iz=1:nz
     cosa             = cos(im(j)*tarr(it)-in(j)*zarr(iz));
     sina             = sin(im(j)*tarr(it)-in(j)*zarr(iz));
     
     Rarr{1}(:,it,iz) = Rarr{1}(:,it,iz) + (Rac(j) + fac{j}{1}.*(Rbc(j)-Rac(j)) )*cosa;
     
     Rarr{2}(:,it,iz) = Rarr{2}(:,it,iz) + fac{j}{2}.*(Rbc(j)-Rac(j))*cosa;
     
     Rarr{3}(:,it,iz) = Rarr{3}(:,it,iz) - im(j)*(Rac(j) + fac{j}{1}.*(Rbc(j)-Rac(j)) )*sina;
    end
  end
end


jacobian = Rarr{2}*rpol*rtor;  %this is (dR/ds)*rpol*rtor

