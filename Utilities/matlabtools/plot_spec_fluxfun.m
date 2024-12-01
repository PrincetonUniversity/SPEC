function plot_spec_fluxfun(data,lvol,ns,nt,z0,ncont,newfig)

%
% PLOT_SPEC_OUTFLUXFUN( DATA, NS, NT, Z0, NCONT, NEWFIG )
% =======================================================
%
% Plots iso-contours of Az on a given cross-section
%
% INPUT
% -----
%   -data     : must be produced by calling read_spec(filename)
%   -lvol     : volume number
%   -ns       : radial resolution for construction of Az
%   -nt       : poloidal resolution for construction of Az
%   -z0       : toroidal angle at which Az is evaluated
%   -ncont    : number of iso-contour lines
%   -newfig   : flag for whether a new figure should be open (=1) or not(=0)
%
% Note: should only work for free-boundary equilibria
% Note: poloidal flux function is psi=-A_phi if A_phi is the covariant component and psi=-R*A_phi if A_phi is the canonical component
%
% written by J.Loizu (2018)


    % Check input
    if lvol>data.output.Mvol || lvol<1
        error('Invalid lvol')
    end

    if ns<1
        error('Invalid ns')
    end

    if nt<1
        error('Invalid nt')
    end

    if ncont<1
        error('Invalid ncont')
    end

    switch newfig
        case 0
            hold on
        case 1
            figure
            hold on
        case 2
            hold off
        otherwise
            error('Invalid newfig')
    end

    % Generate coordinate arrays
    sarr = linspace(-1,1,ns);
    tarr = linspace(0,2*pi,nt);

    acov = get_spec_vecpot(data,lvol,sarr,tarr,z0);
    R = get_spec_R_derivatives(data,lvol,sarr,tarr,z0,'R');
    Z = get_spec_R_derivatives(data,lvol,sarr,tarr,z0,'Z');

    ffun = -acov{2}; 

    if(newfig==1)
     figure; hold on; 
    end

    contour(R{1},Z{1},ffun,ncont,'k')
    
end
