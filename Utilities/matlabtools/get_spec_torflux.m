function psitor = get_spec_torflux(data,lvol,zeta,start,send,ns,nt)

%
% GET_SPEC_TORFLUX( DATA, LVOL, ZETA, START, SEND, NS, NT )
% =========================================================
%
% Computes total enclosed toroidal flux in the poloidal cross-section defined by zeta,
% inside the volume number lvol and across the radial extension defined by start and send
%
% INPUT
% -----
%   -data   : must be produced by calling read_spec(filename)
%   -lvol    : volume number
%   -zeta    : toroidal angle at which the flux is calculated
%   -start   : first point in the radial direction (e.g. start=-1)
%   -send    : last point in the radial direction (e.g. send=+1)
%   -ns      : radial resolution   (e.g. 64)
%   -nt      : poloidal resolution (e.g. 64)
%
% OUPUT
% -----
%   -psitor  : total enclosed toroidal flux
%
%   written by J.Loizu (2016)
%   modified by J.Loizu (01.2017)
%   modified by J.Loizu (06.2017)
%   modified by A.Baillod (06.2019) - added switch for geometry

    %Check inputs
    Igeometry = data.input.physics.Igeometry;
    if lvol==1 && Igeometry~=1 && start==1
        error('InputError: start should be >1 in first volume')
    end

    Mvol = data.output.Mvol;
    if lvol<1 || lvol>Mvol
        error('InputError: Invalid lvol')
    end

    if start<-1 || start>send
        error('InputError: invalid start')
    end

    if send<start || send>1
        error('InputError: invalid send')
    end

    if ns<1
        error('InputError: invalid ns')
    end

    if nt<1
        error('InputError: invalid nt')
    end

    % Prepare coordinate arrays
    sarr = linspace(start,send,ns);
    tarr = linspace(0,2*pi,nt);

    ds   = sarr(2)-sarr(1);

    dth  = tarr(2)-tarr(1);

    if(ds==0 || dth==0)

        psitor = 0;

    else

        Bcontrav = get_spec_magfield(data,lvol,sarr,tarr,zeta);
        jac      = get_spec_jacobian(data,lvol,sarr,tarr,zeta);


        % Compute surface integral

        Bzeta    = Bcontrav{3};
        psitor   = trapz( sarr, trapz( tarr, jac.*Bzeta, 2 ) );
    end

end
