function T = get_spec_polynomial_basis(data,lvol,sarr)

%
% GET_SPEC_POLYNOMIAL_BASIS( DATA, LVOL, SARR )
% =============================================
%
% Return the Chebychev of the Zernike polynomial basis
%
% INPUT
% -----
%   data: generated via read_spec_filename
%   lvol: volume number
%   sarr: radial coordinate, shape (ns, 1)
%
% OUTPUT
% ------
%   T{i}{j}: Polynom of order i (j=1) and its derivative (j=2)
%
%
% Written by A. Baillod (2020)
%

Lrad = double(data.input.physics.Lrad(lvol));
ns   = length(sarr);
Mpol = data.input.physics.Mpol;

Lsingularity = false;
if ((lvol==1) && (data.input.physics.Igeometry~=1))
  Lsingularity = true;
end


T       = cell(Lrad+1,2);  % allocate data for Chebyshev polynomials and their derivatives

T{1}{1} = ones(ns,1);
T{1}{2} = zeros(ns,1);

T{2}{1} = sarr;
T{2}{2} = ones(ns,1);


if( Lsingularity ) % Build zernike polynomials
% Copy pasted from Frotran source, translated to Matlab language
  zernike = zeros(Lrad+1,Mpol+1,3,length(sarr));
  rm  = zeros(size(sarr));
  rm1 = zeros(size(sarr));

  sbar = (1 + sarr ) / 2.0;
  ns = length(sbar);

  for m = 0:Mpol
    if (Lrad >= m)
      zernike(m+1,m+1,2,:) = rm;
      zernike(m+1,m+1,3,:) = double(m)*rm1;
    end

    if (Lrad >= m+2)
      zernike(m+3,m+1,2,:) = double(m+2)     *rm.*sbar.^2 - double(m+1)    *rm;
      zernike(m+3,m+1,3,:) = double((m+2)^2) *rm.*sbar    - double((m+1)*m)*rm1;
    end

    for n = m+4:2:Lrad
      factor1 = double(n) / double(n^2 - m^2);
      factor2 = double(4 * (n-1));
      factor3 = double((n-2+m)^2)/double(n-2) + double((n-m)^2)/double(n);
      factor4 = double((n-2)^2-m^2) / double(n-2);
 
      zernike(n+1, m+1, 2, :) = factor1 * ((factor2*sbar.^2 - factor3) .* reshape(zernike(n-1, m+1, 2, :), [ns,1]) - factor4 * reshape(zernike(n-3, m+1, 2, :), [ns, 1]));
      zernike(n+1, m+1, 3, :) = factor1 * (2.0*factor2*sbar .* reshape(zernike(n-1, m+1, 2, :), [ns,1]) + (factor2*sbar.^2 - factor3) .* reshape(zernike(n-1, m+1, 3, :), [ns,1]) - factor4 * reshape(zernike(n-3, m+1, 3, :), [ns,1]));
    end
 
    rm1 = rm;
    rm  = rm .* sbar;
  end

  for n = 2:2:Lrad
    zernike(n+1,1,1,:) = zernike(n+1,1,1, :) - (-1)^(n/2);
  end

  if (Mpol >= 1)
    for n = 3:2:Lrad
      zernike(n+1,2,1,:) = reshape(zernike(n+1,2,1,:), [ns,1]) - (-1)^((n-1)/2) * double((n+1)/2) *     sbar;
      zernike(n+1,2,2,:) = reshape(zernike(n+1,2,2,:), [ns,1]) - (-1)^((n-1)/2) * double((n+1)/2);
    end
  end

  for m = 0:Mpol
    for n = m:2:Lrad
      zernike(n+1,m+1,:,:) = zernike(n+1,m+1,:,:) / double(n+1);
    end
  end

  % Store in T{}{} structure
  for ll = 1:Lrad+1
    T{ll}{1} = zernike(ll,:,1,:);
    T{ll}{2} = zernike(ll,:,2,:) / 2.0;
  end

else % Otherwise construct Chebychev basis

  for l=3:Lrad+1
    T{l}{1} = 2*sarr.*T{l-1}{1} - T{l-2}{1};
    T{l}{2} = 2*T{l-1}{1} + 2*sarr.*T{l-1}{2} - T{l-2}{2};
  end
end
