function plot_spec_poincare_3d( data, varargin )
%
% PLOT_SPEC_POINCARE_3D( DATA )
% -----------------------------
%
% Plot parts of the 3d stellarator geometry, the Poincare plot and traces a
% field line 
%
% INPUTS
% ------
%  * data:     Obtained with read_spec(filename)
%  * varargin: optional inputs. Can be any pairs of
%        'nt': Number of poloidal points. Default is 64
%        'nz': Number of toroidal points. Default is 64
%        'phiend': Plots 3d shape from phi=0 to phiend. Default is 2pi/Nfp
%        'tstart': Start field line tracing at this value of theta. Can be
%                  an array. Default is 0.


    % Load the colors and some data
    col = EPFL_colors;
    Nvol = double(data.input.physics.Nvol);
    Nfp = double(data.input.physics.Nfp);

    % Set up options
    opt.nt = 64;
    opt.nz = 64;
    opt.phiend = 2*pi / max(Nfp,2);
    opt.tstart = 0;

    l = length(varargin);
    if mod(l,2)~=0
        error('Invalid number of arguments')
    end

    for ii=1:l/2
       opt.(varargin{2*ii-1}) = varargin{2*ii}; 
    end

    opt.tstart = reshape(opt.tstart, length(opt.tstart), 1);

    % Define coordinate arrays in real space
    sarr = 1;                       % We plot on the outermost surface, i.e. the plasma boundary
    tarr = linspace(0,2*pi       ,opt.nt);
    zarr = linspace(0,opt.phiend ,opt.nz);

    R = get_spec_R_derivatives(data,Nvol,sarr,tarr,zarr,'R');
    Z = get_spec_R_derivatives(data,Nvol,sarr,tarr,zarr,'Z');

    R = squeeze(R{1});   
    Z = squeeze(Z{1});


    % Construct cartesian corrdinates 

    X = zeros(opt.nt,opt.nz);
    Y = zeros(opt.nt,opt.nz);

    for it=1:opt.nt
     for iz=1:opt.nz
      X(it,iz) = R(it,iz)*cos(zarr(iz));
      Y(it,iz) = R(it,iz)*sin(zarr(iz));
     end
    end


    % Plot

    figure('Position', [200 200 900 700], 'Color', 'w')

    s = size(X);
    c = zeros(s(1), s(2), 3);
    c(:,:,1) = 1;
    mesh(X,Y,Z,c);
    hold on

    axis equal
    shading interp

    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    set(gca, 'FontSize', 18 )

    % Add white surface on plane y=0
    M = makehgtform('xrotate',pi/2);
    t = hgtransform('Matrix', M);

    pshape = polyshape(X(:,1), Z(:,1));
    plot(pshape, 'Parent', t, 'FaceColor', 'w', 'FaceAlpha', 0.9 )




    % Add Poincare plot
    X = squeeze(data.poincare.R(:,1,:));

    s = size(X);

    Y = zeros(size(X));
    Z = squeeze(data.poincare.Z(:,1,:));
    for i=1:s(1)     %for each field line trajectory
        scatter3(X(i,:),Y(i,:),Z(i,:),10,'.k')
        hold on;
    end

    % Add KAM surfaces
    mn              = data.output.mn;
    im              = data.output.im;
    Rbcmn           = data.output.Rbc;
    Rbsmn           = data.output.Rbs;
    Zbcmn           = data.output.Zbc;
    Zbsmn           = data.output.Zbs;

    nt = 1024;
    tarr = linspace(0, 2*pi, nt);

    X      = zeros(Nvol,nt);
    Y      = zeros(Nvol,nt);
    for i=1:Nvol
        for k=1:mn
            alpha  = double(im(k))*tarr; %phi = 0
            X(i,:) = X(i,:) + Rbcmn(k,i+1)*cos(alpha) + Rbsmn(k,i+1)*sin(alpha);
            Y(i,:) = Y(i,:) + Zbsmn(k,i+1)*sin(alpha) + Zbcmn(k,i+1)*cos(alpha);
        end
    end

    for i=1:size(X,1)
     scatter3(X(i,:),zeros(nt,1),Y(i,:),3,'filled','MarkerFaceColor','r','MarkerEdgeColor','r')
     hold on
    end

    % Add a field line
    tarr = opt.tstart;


    if ~isempty(tarr)
        phi = 0;

        nstep = 1024;
        dstep = opt.phiend/1024;

        for istep=1:nstep
           phi(istep+1) = phi(istep) + dstep;

           B = get_spec_magfield( data, Nvol, 1, tarr(:,istep), phi(istep) );

           Bt = reshape(B{2}, length(B{2}), 1);
           Bz = reshape(B{3}, length(B{2}), 1);

           tarr(:,istep+1) = tarr(:,istep) + Bt./Bz * dstep;

        end

        nline = length(opt.tstart);
        for jj=1:nline
            X = [];
            Y = [];
            Z = [];
            for ii = 1:nstep
                Rd = get_spec_R_derivatives(data,Nvol,sarr,tarr(jj,ii),phi(ii),'R');
                Zd = get_spec_R_derivatives(data,Nvol,sarr,tarr(jj,ii),phi(ii),'Z');

                X(ii) = squeeze(Rd{1}) * cos(phi(ii));
                Y(ii) = squeeze(Rd{1}) * sin(phi(ii));
                Z(ii) = squeeze(Zd{1});
            end

            scatter3( X, Y, Z, 10, 'MarkerFaceColor', col.Leman, 'MarkerEdgeColor', col.Leman)
        end
    end


end