function plot_spec_kam(filename,zetaov2pi,newfig)

% Produces a "Poincare plot" of the KAM surfaces in toroidal geometry
%
% INPUT
%   -filename   : SPEC output hdf5 file
%   -zetaov2pi  : shows the toroidal plane at zeta=2*pi*(zetaov2pi)
%   -newfig     : opens(=1) or not(=0) a new figure
%
%   written by J.Loizu (2016)
%   upgraded by J.Loiyu (07.2017)


Nvol   = h5read(filename,'/Nvol');
tflux  = h5read(filename,'/tflux');
mn     = h5read(filename,'/mn');
im     = h5read(filename,'/im');
in     = h5read(filename,'/in');
Rbcmn  = h5read(filename,'/Rbc');
Rbsmn  = h5read(filename,'/Rbs');
Zbcmn  = h5read(filename,'/Zbc');
Zbsmn  = h5read(filename,'/Zbs');


% Compute (x,y) coordinates of each KAM surface

zeta  = zetaov2pi*(2*pi);

nth    = 2048;
dth    = 2*pi/nth;
theta  = dth:dth:2*pi; 

X      = zeros(Nvol,nth);
Y      = zeros(Nvol,nth);

for i=1:Nvol
 for k=1:mn
 alpha  = double(im(k))*theta-double(in(k))*zeta;
 X(i,:) = X(i,:) + Rbcmn(k,i+1)*cos(alpha) + Rbsmn(k,i+1)*sin(alpha);
 Y(i,:) = Y(i,:) + Zbsmn(k,i+1)*sin(alpha) + Zbcmn(k,i+1)*cos(alpha);
 end
end


% Plot Poincare section

if(newfig==1)
figure
end
hold on

for i=1:Nvol
 scatter(X(i,:),Y(i,:),3,'filled','r')
end
axis equal
hold on
set(gca,'FontSize',12)
xlabel('R','FontSize',12)
ylabel('Z','FontSize',12)
%xlim([-1.1*rmax 1.1*rmax])
%ylim([-1.1*zmax 1.1*zmax])
