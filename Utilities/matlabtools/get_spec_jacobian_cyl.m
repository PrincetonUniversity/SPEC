function jacobian = get_spec_jacobian_cyl(data,lvol,sarr,tarr,zarr)
 
% Calculates Jacobian in volume number lvol in cylindrical coordinates
%
% INPUT
%   -data      : must be produced by calling e.g. read_spec_grid(filename)
%   -lvol      : volume number
%   -sarr      : is the array of values for the s-coordinate
%   -tarr      : is the array of values for the theta-coordinate
%   -zarr      : is the array of values for the zeta-coordinate
%
% OUTPUT
%   -jacobian  : Jacobian of the coordinates with size ns*nt*nz where 
%                ns=length(sarr),nt=length(zarr),nt=length(zarr)
%
% Note: Stellarator symmetry is assumed
%
% written by A.Baillod (2019)


Rac     = data.Rbc(:,lvol);   % inner volume boundary harmonics
Rbc     = data.Rbc(:,lvol+1); % outer volume boundary harmonics

sarr    = transpose(sarr);
ns      = length(sarr);
nt      = length(tarr);
nz      = length(zarr);
sbar    = (sarr+1)/2;

mn      = data.mn;
im      = double(data.im);
in      = double(data.in);

Rarr    = cell(2); % allocate data for R-array and its derivatives in (s,theta)

for k=1:2
  Rarr{k} = zeros(ns,nt,nz);
end

fac     = cell(mn,2);      % allocate data for regularization factors and their derivatives

% Construct regularization (for lvol=1) factors and their derivatives

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
        fac{j}{2}  = 1.0 / 2.0;
    end
    
    
    
%     
%   if(lvol>1 || im(j)==0) 
%    fac{j}{1}  = sqrt(sbar);
%    fac{j}{2}  = 0.25./sqrt(sbar);
%   else
%    fac{j}{1}  = sbar.^(im(j)/2);
%    fac{j}{2}  = (im(j)/4)*sbar.^(im(j)/2-1);
%  end
end



% Construct (R,Z) arrays and their derivatives

for j=1:mn
  for it=1:nt
    for iz=1:nz
     cosa = cos(im(j)*tarr(it)-in(j)*zarr(iz));
     
     Rarr{1}(:,it,iz) = Rarr{1}(:,it,iz) + (Rac(j) + fac{j}{1}.*(Rbc(j)-Rac(j)) )*cosa;
     
     Rarr{2}(:,it,iz) = Rarr{2}(:,it,iz) +           fac{j}{2}.*(Rbc(j)-Rac(j))  *cosa;
    end
  end
end


jacobian = Rarr{1}.*Rarr{2};
