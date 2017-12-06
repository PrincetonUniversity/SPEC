function gmat = get_spec_metric(data,lvol,sarr,tarr,zarr)
 
 
% Calculates covariant metric elements of the coordinates in volume number lvol
%
% INPUT
%   -data   : must be produced by calling e.g. read_spec_grid(filename)
%   -lvol   : volume number
%   -sarr   : is the array of values for the s-coordinate
%   -tarr   : is the array of values for the theta-coordinate
%   -zarr   : is the array of values for the zeta-coordinate
%
% OUTPUT
%   -gmat   : 3x3 matrix of metric coefficients 
%
% Note: Stellarator symmetry is assumed
%
% written by J.Loizu (2016)


Rac     = data.Rbc(:,lvol);   % inner volume boundary harmonics
Zas     = data.Zbs(:,lvol);
Rbc     = data.Rbc(:,lvol+1); % outer volume boundary harmonics
Zbs     = data.Zbs(:,lvol+1);

sarr    = transpose(sarr);
ns      = length(sarr);
nt      = length(tarr);
nz      = length(zarr);
sbar    = (sarr+1)/2;

mn      = data.mn;
im      = double(data.im);
in      = double(data.in);

Rarr    = cell(4); % allocate data for R-array and its derivatives in (s,theta,zeta)
Zarr    = cell(4); % allocate data for Z-array and its derivatives in (s,theta,zeta)

for k=1:4
Rarr{k} = zeros(ns,nt,nz); 
Zarr{k} = zeros(ns,nt,nz); 
end

fac     = cell(mn,2); % allocate data for regularization factors and their derivatives

gmat    = cell(3,3);  % allocate data for the metric matrix

for k=1:3
 for p=1:3
  gmat{k}{p} = zeros(ns,nt,nz); 
 end
end


% Construct regularization (for lvol=1) factors and their derivatives

for j=1:mn
  if(lvol>1 || im(j)==0) 
   fac{j}{1}  = sbar;
   fac{j}{2}  = 0.5*ones(ns,1);
  else
   fac{j}{1}  = sbar.^(im(j)/2);
   fac{j}{2}  = (im(j)/4)*sbar.^(im(j)/2-1);
  end
end


% Construct (R,Z) arrays and their derivatives

for j=1:mn
  for it=1:nt
    for iz=1:nz
     cosa = cos(im(j)*tarr(it)-in(j)*zarr(iz));
     sina = sin(im(j)*tarr(it)-in(j)*zarr(iz));
     Rarr{1}(:,it,iz) = Rarr{1}(:,it,iz) + (Rac(j) + fac{j}{1}.*(Rbc(j)-Rac(j)) )*cosa;
     Zarr{1}(:,it,iz) = Zarr{1}(:,it,iz) + (Zas(j) + fac{j}{1}.*(Zbs(j)-Zas(j)) )*sina;
     
     Rarr{2}(:,it,iz) = Rarr{2}(:,it,iz) + fac{j}{2}.*(Rbc(j)-Rac(j))*cosa;
     Zarr{2}(:,it,iz) = Zarr{2}(:,it,iz) + fac{j}{2}.*(Zbs(j)-Zas(j))*sina;
     
     Rarr{3}(:,it,iz) = Rarr{3}(:,it,iz) - im(j)*(Rac(j) + fac{j}{1}.*(Rbc(j)-Rac(j)) )*sina;
     Zarr{3}(:,it,iz) = Zarr{3}(:,it,iz) + im(j)*(Zas(j) + fac{j}{1}.*(Zbs(j)-Zas(j)) )*cosa;
     
     Rarr{4}(:,it,iz) = Rarr{4}(:,it,iz) + in(j)*(Rac(j) + fac{j}{1}.*(Rbc(j)-Rac(j)) )*sina;
     Zarr{4}(:,it,iz) = Zarr{4}(:,it,iz) - in(j)*(Zas(j) + fac{j}{1}.*(Zbs(j)-Zas(j)) )*cosa;
    end
  end
end


% Construct metric elements

gmat{1}{1} = Rarr{2}.^2 + Zarr{2}.^2;                %gss
gmat{2}{2} = Rarr{3}.^2 + Zarr{3}.^2;                %gtt
gmat{3}{3} = Rarr{1}.^2 + Rarr{4}.^2 + Zarr{4}.^2;   %gzz
gmat{1}{2} = Rarr{2}.*Rarr{3} + Zarr{2}.*Zarr{3};    %gst
gmat{1}{3} = Rarr{2}.*Rarr{4} + Zarr{2}.*Zarr{4};    %gsz
gmat{2}{3} = Rarr{3}.*Rarr{4} + Zarr{3}.*Zarr{4};    %gtz

gmat{2}{1} = gmat{1}{2};  % by symmetry of g
gmat{3}{1} = gmat{1}{3};
gmat{3}{2} = gmat{2}{3};


