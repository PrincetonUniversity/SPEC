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
ns = length(sarr);

T       = cell(Lrad+1,2);  % allocate data for Chebyshev polynomials and their derivatives

T{1}{1} = ones(ns,1);
T{1}{2} = zeros(ns,1);

T{2}{1} = sarr;
T{2}{2} = ones(ns,1);


for l=3:Lrad+1
  T{l}{1} = 2*sarr.*T{l-1}{1} - T{l-2}{1};
  T{l}{2} = 2*T{l-1}{1} + 2*sarr.*T{l-1}{2} - T{l-2}{2};
end

