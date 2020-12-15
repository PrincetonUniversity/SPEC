function fac = get_spec_regularization_factor(data, lvol, sarr, ForG)
%
%
% Computes the regularisation factor in the correct geometry
%
% INPUT
% -----
%   data:     Produced by fdata_from_data(data);
%   lvol:     Volume number
%   sarr:     s-coordinate array, shape (ns, 1)
%   ForG:     'F' for field reg. factor, 'G' for geometry reg. factor
%
% OUPUT
% -----
%   fac: mnx2 cell array composed of fj and its derivatives
%
%
% Written by A.Baillod (2019)
%
%

sarr = transpose(sarr);

Igeometry= data.input.physics.Igeometry;
mn      = data.output.mn;
im      = double(data.output.im);
ns      = length(sarr);
fac     = cell(mn,2);
sbar    = (1+sarr)/2.0;
Mregular= double(data.input.numerics.Mregular);


regumm = im / 2.0;
if Mregular>1
   ind = find(regumm>Mregular);
   regumm(ind) = Mregular / 2.0;
end

halfmm = im / 2.0;



if ForG=='G'
    switch Igeometry
        case 1 % Slab geometry

            for j=1:mn
                fac{j}{1}  = sbar;
                fac{j}{2}  = 0.5*ones(ns,1);
            end

        case 2 % Cylindrical geometry
            %TODO: add glab if Lextrap is not set to 1!
            for j=1:mn
                if(lvol==1) 
                    if im(j)==0
                        fac{j}{1} = sqrt(sbar);
                        fac{j}{2} = 0.25 * fac{j}{1} ./ sbar;
                    else
                        fac{j}{1} = sbar.^halfmm(j); 
                        fac{j}{2} = 0.5 * halfmm(j) * fac{j}{1} ./ sbar;
                    end
                else
                    fac{j}{1}  = sbar;
                    fac{j}{2}  = 0.5 * ones(ns,1);
                end
            end

        case 3 % Toroidal geometry

            for j=1:mn
                if lvol==1 %coordinate singularity
                    if im(j)==0
                       fac{j}{1} = sbar; 
                       fac{j}{2} = 0.5;
                    else
                       fac{j}{1} = sbar.^halfmm(j);
                       fac{j}{2} = 0.5 * halfmm(j) * fac{j}{1} ./ sbar;
                    end
                else
                    fac{j}{1}  = sbar;
                    fac{j}{2}  = 0.5 * ones(ns,1);
                end
            end
        otherwise
            error('Unsupported geometry in get_spec_regularisation_factor')
    end
    
    
    
elseif ForG=='F'
    switch Igeometry
        case 1
            for j=1:mn
                fac{j}{1}  = ones(ns,1);
                fac{j}{2}  = zeros(ns,1);
            end
            
        case 2
            for j=1:mn
                if lvol==1 % coordinate singularity
                   if im(j)==0
                      fac{j}{1} = 1;
                      fac{j}{2} = 0;
                   else
                      fac{j}{1} = sbar.^regumm(j);
                      fac{j}{2} = 0.5 * regumm(j) * fac{j}{1} ./ sbar;
                   end
                else
                  fac{j}{1} = 1;
                  fac{j}{2} = 0;
                end
            end
            
        case 3
            for j=1:mn
                if lvol==1 % coordinate singularity
                   if im(j)==0
                      fac{j}{1} = 1;
                      fac{j}{2} = 0;
                   else
                      fac{j}{1} = sbar.^regumm(j);
                      fac{j}{2} = 0.5 * regumm(j) * fac{j}{1} ./ sbar;
                   end
                else
                  fac{j}{1} = 1;
                  fac{j}{2} = 0;
                end
            end   
          
        otherwise
            error('Unsupported geometry in get_spec_regularisation_factor')
    end
    
else
    error('Unsupported ForG value in get_spec_regularisation_factor.')
end
