function plot_spec_kam(data,zetaov2pi,newfig)

% Produces a "Poincare plot" of the KAM surfaces.
%
% INPUT
%   -data       : obtained from read_spec(fname)
%   -zetaov2pi  : shows the toroidal plane at zeta=2*pi*(zetaov2pi)
%   -newfig     : opens(=1) or not(=0) a new figure, or overplots(=2) on existing figure
%
%   written by J.Loizu (2016)
%   upgraded by J.Loizu (07.2017)
%   modified by A. Baillod (06.2019)
%   modified by J.Loizu (01.2020)


Nvol            = double(data.input.physics.Nvol);
mn              = data.output.mn;
im              = data.output.im;
in              = data.output.in;
Rbcmn           = data.output.Rbc;
Rbsmn           = data.output.Rbs;
Zbcmn           = data.output.Zbc;
Zbsmn           = data.output.Zbs;
Igeometry       = data.input.physics.Igeometry;

% Compute (x,y) coordinates of each KAM surface

zeta  = zetaov2pi*(2*pi);

nth    = 2048;
dth    = 2*pi/nth;
theta  = dth:dth:2*pi; 

X      = zeros(Nvol,nth);
Y      = zeros(Nvol,nth);


switch Igeometry
    case 1
        for i=1:Nvol
            X(i,:) = theta;
            for k=1:mn
                alpha  = double(im(k))*theta-double(in(k))*zeta;
                Y(i,:) = Y(i,:) + Rbcmn(k,i+1)*cos(alpha) + Rbsmn(k,i+1)*sin(alpha);
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

for i=1:Nvol
 scatter(X(i,:),Y(i,:),3,'filled','r')
 hold on
end

if Igeometry~=1
    axis equal
end
hold on
set(gca,'FontSize',12)
xlabel('R','FontSize',12)
ylabel('Z','FontSize',12)
