function jacobian = get_spec_jacobian(data,lvol,sarr,tarr,zarr)
 
 
% Calculates the coordinate Jacobian in volume number lvol
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
% written by J.Loizu (2016)
% modified by A. Baillod (2019)

[sarr, Rarr] = get_spec_R_derivatives(data, lvol, sarr, tarr, zarr, 'R');
[sarr, Zarr] = get_spec_R_derivatives(data, lvol, sarr, tarr, zarr, 'Z');


switch data.Igeometry
    case 1
        jacobian = Rarr{2};
    case 2
        jacobian = Rarr{1}.*Rarr{2};
    case 3
        jacobian = Rarr{1}.*(Rarr{3}.*Zarr{2} - Rarr{2}.*Zarr{3});
    otherwise
        error('Unsupported geometry in get_spec_jacobian')
end

