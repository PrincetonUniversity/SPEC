function jacobian = get_spec_jacobian(data,lvol,sarr,tarr,zarr)
 
%
% GET_SPEC_JACOBIAN( DATA, LVOL, SARR, TARR, ZARR )
% =================================================
%
% Calculates the coordinate Jacobian in volume number lvol
%
% INPUT
% -----
%   -data      : must be produced by calling read_spec(filename)
%   -lvol      : volume number
%   -sarr      : is the array of values for the s-coordinate
%   -tarr      : is the array of values for the theta-coordinate
%   -zarr      : is the array of values for the zeta-coordinate
%
% OUTPUT
% ------
%   -jacobian  : Jacobian of the coordinates with size ns*nt*nz where ns=length(sarr),nt=length(zarr),nt=length(zarr)
%
% Note: Stellarator symmetry is assumed
%
% written by J.Loizu (2016)
% modified by A. Baillod (2019)

    % Check input
    Istellsym = data.input.physics.Istellsym;
    if Istellsym==0
       error('Non stellarator symmetric not implemented') 
    end

    Rarr = get_spec_R_derivatives(data, lvol, sarr, tarr, zarr, 'R');

    switch data.input.physics.Igeometry
        case 1
            jacobian = Rarr{2};
        case 2
            jacobian = Rarr{1}.*Rarr{2};
        case 3
            Zarr = get_spec_R_derivatives(data, lvol, sarr, tarr, zarr, 'Z');
            jacobian = Rarr{1}.*(Rarr{3}.*Zarr{2} - Rarr{2}.*Zarr{3});
        otherwise
            error('Unsupported geometry in get_spec_jacobian')
    end

end