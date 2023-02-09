%% extract_shear(DATA, n, shift, n_surf)
% =======================================
%
% Compute the shear (from iota) out of 
% a SPEC out data
%
% INPUT
% -----
%   -data :   must be produced using read_spec
%   -n :      length(iota) - (n_surf+shift-1);
%   -shift :  shift of the indices we take into account when computing
%             shear (1st indices might not be relevant (too close to axis))
%   -n_surf : index of the first surface that encircles the magnetic axis
%
% ------------------------------------%
% Written by S.Guinchard (05/12/22)   % 
% ------------------------------------%
function out = extract_shear(d, n, shift, n_surf)
    
  id = n_surf+shift;
  radial_coord     = (d.transform.fiota(:,1));       % extract radial coordinate
  out.mat_iota     =  d.transform.fiota(id:end,2);   % extract iota
  out.mat_r_coord  = radial_coord(id:end);           % truncates the radial coordinates (remove 5 first terms)
  out.mat_s_coord  = ((out.mat_r_coord + 1)./2);     % change of variable r <--> s
  out.scan_11      = d.output.Rbc(11,2);             % Value of R11
  out.scan_10      = d.output.Rbc(2,2);              % Value of R10
 
  % CENTERED FINITE DIFFERENCES 
    out.derivatives(1) = (out.mat_iota(2) - out.mat_iota(1))   / (out.mat_s_coord(2) - out.mat_s_coord(1));
    out.derivatives(n) = (out.mat_iota(n) - out.mat_iota(n-1)) / (out.mat_s_coord(n) - out.mat_s_coord(n-1));
    
 for j=2:n-1
    out.derivatives(j) = (out.mat_iota(j+1) - out.mat_iota(j-1)) / (out.mat_s_coord(j+1) - out.mat_s_coord(j-1));
 end
 
 out.coeff      =    out.mat_s_coord ./ (out.mat_iota) ;
 out.shear      =    (out.coeff)'     .* out.derivatives;
 out.avg_shear  =    mean(out.shear);                     % Avg shear 
end