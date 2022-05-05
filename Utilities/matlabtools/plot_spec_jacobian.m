function rzbdata = plot_spec_jacobian(data,lvol,sarr,tarr,zarr,newfig)

%
% PLOT_SPEC_JACOBIAN( DATA, LVOL, SARR, TARR, ZARR, NEWFIG )
% ==========================================================
%
% Produces plot of sqrt(g) in (R,Z,zarr) cross-section(s)
%
% INPUT
% -----
%   -data    : data obtained via read_spec(filename)
%   -lvol    : volume number. Set to 0 for plotting all volumes
%   -sarr    : is the array of values for the s-coordinate ('d' for default)
%   -tarr    : is the array of values for the theta-coordinate ('d' for default)
%   -zarr    : is the array of values for the zeta-coordinate ('d' for default)
%   -newfig  : opens(=1) or not(=0) a new figure, or overwrite existing one
%   (=2)
%
% OUTPUT
% ------
%   -rzbdata : cell structure with 3 arrays: R-data, Z-data, |B|-data
%
% written by J.Loizu (2016)

if(sarr=='d')
sarr=linspace(-1,1,64);
end

if(tarr=='d')
tarr=linspace(0,2*pi,64);
end

if(zarr=='d')
zarr=0;
end

% Check input
if (length(sarr)>1) && length(tarr)>1 && length(zarr)>1
   error('This is a 2d plotting routine; one input array has to be a scalar') 
end

switch newfig
    case 0
        hold on
    case 1
        figure
        hold on
    case 2
        hold off
end

% Allocate memory
rzbdata = cell(3);

if lvol==0
    lstart=1;
    lend  =data.output.Mvol;
else
    lstart=lvol;
    lend  =lvol;
end

for ivol=lstart:lend
    % Compute sqrt(g)
    jac   = squeeze(get_spec_jacobian(data,ivol,sarr,tarr,zarr));

    % Compute function (R,Z)(s,theta,zeta)
    R = get_spec_R_derivatives(data,ivol,sarr,tarr,zarr,'R');
    Z = get_spec_R_derivatives(data,ivol,sarr,tarr,zarr,'Z');

    R = R{1};   
    Z = Z{1};

    % Plot
    Rtemp = R;
    Ztemp = Z;
    switch data.input.physics.Igeometry
        case 1
            R = tarr;
            Z = Rtemp;
        case 2
            for it=1:length(tarr)
                R(:,it,:) = Rtemp(:,it,:) .* cos(tarr(it));
                Z(:,it,:) = Rtemp(:,it,:) .* sin(tarr(it));
            end
        case 3
            R = squeeze(Rtemp);
            Z = squeeze(Ztemp);
    end


    pcolor(R,Z,jac); 
    shading interp; 
    colorbar
    hold on
end

axis equal
title('|B|');
xlabel('R');
ylabel('Z');

% Output data

rzbdata{1} = R;
rzbdata{2} = Z;
rzbdata{3} = jac;
