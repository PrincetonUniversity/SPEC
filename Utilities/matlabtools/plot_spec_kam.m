function plot_spec_kam(pdata,zetaov2pi,newfig)

% Produces a "Poincare plot" of the KAM surfaces in toroidal geometry
% Assumes stellarator symmetry,
%
% INPUT
%   -pdata      : obtained from pdata_from_data(data)
%   -zetaov2pi  : shows the toroidal plane at zeta=2*pi*(zetaov2pi)
%   -newfig     : opens(=1) or not(=0) a new figure
%
%   written by J.Loizu (2016)
%   upgraded by J.Loiyu (07.2017)
%   Modified by A. Baillod (06.2019)


% Nvol   = h5read(filename,'/Nvol');
% mn     = h5read(filename,'/mn');
% im     = h5read(filename,'/im');
% in     = h5read(filename,'/in');
% Rbcmn  = h5read(filename,'/Rbc');
% Rbsmn  = h5read(filename,'/Rbs');
% Zbcmn  = h5read(filename,'/Zbc');
% Zbsmn  = h5read(filename,'/Zbs');
% Igeometry  = h5read(filename,'/Igeometry');

Nvol            = pdata.Mvol - pdata.Lfreebound;
mn              = pdata.mn;
im              = pdata.im;
in              = pdata.in;
Rbcmn           = pdata.Rbc;
Rbsmn           = pdata.Rbs;
Zbcmn           = pdata.Zbc;
Zbsmn           = pdata.Zbs;
Igeometry       = pdata.Igeometry;



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
                X(i,:) = X(i,:) + Rbcmn(k,i+1)*cos(alpha).*cos(theta);
                Y(i,:) = Y(i,:) + Rbcmn(k,i+1)*cos(alpha).*sin(theta);
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
