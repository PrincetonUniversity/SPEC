%% get_metric_QH(BDATA, C)
% ======================
%
% Extract the QH metric from a Booz_xForms output
%
% INPUT
% -----
%   -bdata: must be produced calling read_boozer
%   -c    : constant s.t c*m = n (Fourier modes)
%   
% OUTPUT
% ------
%   -metric: the value of the QH_metric
%
% ------------------------------------%
% Written by S.Guinchard (05/12/22)   % 
% ------------------------------------%
function metric = get_metric_QH(b, c)

bmnc_b = b.Booz_xForms.Outputs.bmnc_b;
xn_b   = b.Booz_xForms.Outputs.xn_b;
xm_b   = b.Booz_xForms.Outputs.xm_b;
nfp    = double(b.Booz_xForms.Inputs.nfp);

Mpol   = b.Booz_xForms.Inputs.mpol;
Ntor   = b.Booz_xForms.Inputs.ntor;
eq     = find(xn_b ~= c*xm_b);       % indices s.t c*m ~= n
m      = find(abs(xm_b) <= Mpol );   % indices s.t abs(m) < Mpol
n      = find(abs(xn_b) <= Ntor*nfp);    % indices s.t abs(n) < Ntor
[val, indm, indn]  = intersect(m,n); % val = indices s.t abs(m) < Mpol & abs(n) < Ntor
ind_to_sum = intersect(val,eq);      % ind_to_sum = indices s.t c*m=!n and previous conditions

metric = sum(bmnc_b(ind_to_sum).^2)/bmnc_b(1)^2; 


end