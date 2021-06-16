function plot_spec_kam(data, nz0, newfig)

%
% PLOT_SPEC_KAM( DATA, NZ0, NEWFIG )
% ========================================
%
% Produces a "Poincare plot" of the KAM surfaces.
%
% INPUT
% -----
%   -data       : obtained from read_spec(fname)
%   -nz0        : show the toroidal plane number nz0
%   -newfig     : opens(=1) or not(=0) a new figure, or overplots(=2) on existing figure
%
%   written by J.Loizu (2016)
%   upgraded by J.Loizu (07.2017)
%   modified by A. Baillod (06.2019)
%   modified by J.Loizu (01.2020)

Np = double(data.input.physics.Nfp);
%Nplan = double(size(data.poincare.R,2));   % # of toroidal planes
%zetaov2pi = (nz0-1) / (Np * Nplan);
zetaov2pi = 0;

Nvol            = double(data.input.physics.Nvol);
mn              = data.output.mn;
im              = data.output.im;
in              = data.output.in;
Rbcmn           = data.output.Rbc;
Rbsmn           = data.output.Rbs;
Zbcmn           = data.output.Zbc;
Zbsmn           = data.output.Zbs;
Igeometry       = data.input.physics.Igeometry;
try
 rpol           = data.input.physics.rpol;
catch
 rpol           = 1;
end

% Compute (x,y) coordinates of each KAM surface

zeta  = zetaov2pi*(2*pi);

nth    = 2048;
dth    = 2*pi/nth;
theta  = dth:dth:2*pi; 

X      = zeros(Nvol,nth);
Y      = zeros(Nvol,nth);


switch Igeometry
    case 1
        X     = zeros(Nvol+1,nth);
        Y     = zeros(Nvol+1,nth);
        for i=1:Nvol+1
            X(i,:) = rpol*theta;
            for k=1:mn
                alpha  = double(im(k))*theta-double(in(k))*zeta;
                Y(i,:) = Y(i,:) + Rbcmn(k,i)*cos(alpha) + Rbsmn(k,i)*sin(alpha);
            end
        end
    case 2
        for i=1:Nvol
            for k=1:mn
                alpha = double(im(k))*theta-double(in(k))*zeta;
                X(i,:) = X(i,:) + (Rbcmn(k,i+1)*cos(alpha) + Rbsmn(k,i+1)*sin(alpha)).*cos(theta);
                Y(i,:) = Y(i,:) + (Rbcmn(k,i+1)*cos(alpha) + Rbsmn(k,i+1)*sin(alpha)).*sin(theta);
            end
        end
    case 3
        for i=1:Nvol
            for k=1:mn
                alpha  = double(im(k))*theta-double(in(k))*zeta;
                X(i,:) = X(i,:) + Rbcmn(k,i+1)*cos(alpha) + Rbsmn(k,i+1)*sin(alpha);
                Y(i,:) = Y(i,:) + Zbsmn(k,i+1)*sin(alpha) + Zbcmn(k,i+1)*cos(alpha);
            end
        end
    otherwise
        error('Unsupported geometry')
end


% Plot Poincare section

switch newfig
    case 0
        hold on
    case 1
        figure
        hold on
    case 2
        hold off
end

for i=1:size(X,1)
 scatter(X(i,:),Y(i,:),3,'filled','r')
 hold on
end

hold on
set(gca,'FontSize',12)

if Igeometry~=1
 axis equal
 xlabel('R','FontSize',12)
 ylabel('Z','FontSize',12)
else
 xlabel('\theta r_{pol}','FontSize',12)
 ylabel('R','FontSize',12)
end
