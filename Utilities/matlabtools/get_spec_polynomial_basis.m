function T = get_spec_polynomial_basis(data,lvol)

%
% GET_SPEC_POLYNOMIAL_BASIS( DATA, LVOL )
% =======================================
%
% Return the Chebychev of the Zernike polynomial basis
%
% INPUT
% -----
%   data: generated via read_spec_filename
%   lvol: volume number
%
% OUTPUT
% ------
%   T{i}{j}: Polynom of order i (j=1) and its derivative (j=2)
%
%
% Written by A. Baillod (2020)
%

T       = cell(Lrad+1,2);  % allocate data for Chebyshev polynomials and their derivatives

T{1}{1} = ones(ns,1);
T{1}{2} = zeros(ns,1);

T{2}{1} = sarr;
T{2}{2} = ones(ns,1);


if(lvol==1) % Then use Zernike basis
  % TODO to complete
else % Then use Chebychev basis
  for l=3:Lrad+1
    T{l}{1} = 2*sarr.*T{l-1}{1} - T{l-2}{1};
    T{l}{2} = 2*T{l-1}{1} + 2*sarr.*T{l-1}{2} - T{l-2}{2};
  end
end







