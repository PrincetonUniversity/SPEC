function gmat = get_spec_metric(data,lvol,sarr,tarr,zarr)
 
%
% GET_SPEC_METRIC( DATA, LVOL, SARR, TARR, ZARR )
% ===============================================
%
% Calculates covariant metric elements of the coordinates in volume number lvol
%
% INPUT
% -----
%   -data   : must be produced by calling e.g. read_spec(filename)
%   -lvol   : volume number
%   -sarr   : is the array of values for the s-coordinate
%   -tarr   : is the array of values for the theta-coordinate
%   -zarr   : is the array of values for the zeta-coordinate
%
% OUTPUT
% ------
%   -gmat   : 3x3 matrix of metric coefficients 
%
% Note: Stellarator symmetry is assumed
%
% written by J.Loizu (2016)
% modified by A. Baillod (2019)]

    % Check input
    Istellsym = data.input.physics.Istellsym;
    if Istellsym==0
       error('Non stellarator symmetric not implemented') 
    end
    
    % Allocate data for the metric matrix
    ns      = length(sarr);
    nt      = length(tarr);
    nz      = length(zarr);

    gmat    = cell(3,3); 

    for k=1:3
     for p=1:3
      gmat{k,p} = zeros(ns,nt,nz); 
     end
    end

    % Get R and its derivatives
    Rarr = get_spec_R_derivatives(data,lvol,sarr,tarr,zarr,'R');
    % Construct metric elements

    switch data.input.physics.Igeometry
        case 1 % Slab geometry
            gmat{1,1} = Rarr{2}.^2;                    
            gmat{2,2} = 1 + Rarr{3}.^2;                
            gmat{3,3} = 1 + Rarr{4}.^2;                
            gmat{1,2} = Rarr{2}.*Rarr{3};              
            gmat{1,3} = Rarr{2}.*Rarr{4};              
            gmat{2,3} = Rarr{3}.*Rarr{4};              
            
        case 2 % Cylindrical geometry
            gmat{1,1} = Rarr{2}.^2;                             %gss
            gmat{2,2} = Rarr{3}.^2 + Rarr{1}.^2;                %gtt
            gmat{3,3} = Rarr{4}.^2 + 1;                         %gzz
            gmat{1,2} = Rarr{2}.*Rarr{3};                       %gst
            gmat{1,3} = Rarr{2}.*Rarr{4};                       %gsz
            gmat{2,3} = Rarr{3}.*Rarr{4};                       %gtz
            
        case 3 %Toroidal geometry
            % Get Z and its derivatives
            Zarr = get_spec_R_derivatives(data,lvol,sarr,tarr,zarr,'Z');
            
            gmat{1,1} = Rarr{2}.^2 + Zarr{2}.^2;                %gss
            gmat{2,2} = Rarr{3}.^2 + Zarr{3}.^2;                %gtt
            gmat{3,3} = Rarr{1}.^2 + Rarr{4}.^2 + Zarr{4}.^2;   %gzz
            gmat{1,2} = Rarr{2}.*Rarr{3} + Zarr{2}.*Zarr{3};    %gst
            gmat{1,3} = Rarr{2}.*Rarr{4} + Zarr{2}.*Zarr{4};    %gsz
            gmat{2,3} = Rarr{3}.*Rarr{4} + Zarr{3}.*Zarr{4};    %gtz
        otherwise
            error('Unsupported geometry in get_spec_metric.m')
    end

    gmat{2,1} = gmat{1,2};  % by symmetry of g
    gmat{3,1} = gmat{1,3};
    gmat{3,2} = gmat{2,3};


end


