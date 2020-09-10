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
%   sarr: radial coordinate
%
% OUTPUT
% ------
%   T{i}{j}: Polynom of order i (j=1) and its derivative (j=2)
%
%
% Written by A. Baillod (2020)
%

Lrad = data.input.physics.Lrad(lvol);
ns   = length(sarr);
Mpol = data.input.physics.Mpol;

Lsingularity = false;
if ((lvol==1) && (data.input.Igeometry~=1))
  Lsingularity = true;
end


T       = cell(Lrad+1,2);  % allocate data for Chebyshev polynomials and their derivatives

T{1}{1} = ones(ns,1);
T{1}{2} = zeros(ns,1);

T{2}{1} = sarr;
T{2}{2} = ones(ns,1);


if( Lsingularity ) % Build zernike polynomials
% Copy pasted from Frotran source, translated to Matlab language
  zernike = zeros(Lrad+1,Mpol+1,2);
  rm = 0.0;
  rm1 = 0.0;

  for m = 0:mpol
    if (Lrad >= m)
      zernike(m,m,0) = rm
      zernike(m,m,1) = double(m)*rm1
    end

    if (lrad >= m+2)
      zernike(m+2,m,0) = double(m+2)     *rm*r^2 - double(m+1)    *rm
      zernike(m+2,m,1) = double((m+2)^2)*rm*r    - double((m+1)*m)*rm1
    end

    for n = m+4:2:Lrad
      factor1 = double(n) / double(n^2 - m^2)
      factor2 = double(4 * (n-1))
      factor3 = double((n-2+m)^2)/double(n-2) + double((n-m)^2)/double(n)
      factor4 = double((n-2)^2-m^2) / double(n-2)
 
      zernike(n, m, 0) = factor1 * ((factor2*r^2 - factor3)*zernike(n-2,m    ,0) - factor4*zernike(n-4,m,0))
      zernike(n, m, 1) = factor1 * (two*factor2*r*zernike(n-2,m,0) + (fact    or2*r^2 - factor3)*zernike(n-2,m,1) - factor4*zernike(n-4,m,1))
    enddo
 
    rm1 = rm
    rm = rm * r
  end

  for n = 2:2:lrad
    zernike(n,0,0) = zernike(n,0,0) - (-1)^(n/2)
  end

  if (mpol >= 1) then
    for n = 3:2:lrad
      zernike(n,1,0) = zernike(n,1,0) - (-1)^((n-1)/2) * double((n+1)/2) *     r
      zernike(n,1,1) = zernike(n,1,1) - (-1)^((n-1)/2) * double((n+1)/2)
    enddo
  end

  for m = 0:mpol
    for n = m:2:lrad
      zernike(n,m,:) = zernike(n,m,:) / double(n+1)
    end
  end

  % Store in T{}{} structure
  for ll = 0:Lrad
    T{ll}{1} = zernike(ll,:,0)      ;
    T{ll}{2} = zernike(ll,:,1) / 2.0;
  end

else % Otherwise construct Chebychev basis

  for l=3:Lrad+1
    T{l}{1} = 2*sarr.*T{l-1}{1} - T{l-2}{1};
    T{l}{2} = 2*T{l-1}{1} + 2*sarr.*T{l-1}{2} - T{l-2}{2};
  end
end
