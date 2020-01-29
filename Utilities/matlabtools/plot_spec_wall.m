function plot_spec_wall(data,zetaov2pi,newfig)

% Produces a "Poincare plot" of the computational boundary surface in toroidal geometry.
%
% INPUT
%   -data       : obtained from read_spec(fname) 
%   -zetaov2pi  : shows the toroidal plane at zeta=2*pi*(zetaov2pi)
%   -newfig     : opens(=1) or not(=0) a new figure
%
% written by J.Loizu (2018)
% modified by J.Loizu (2020)


mn     = data.output.mn;
im     = data.output.im;
in     = data.output.in;
Rbcmn  = data.output.Rbc;
Rbsmn  = data.output.Rbs;
Zbcmn  = data.output.Zbc;
Zbsmn  = data.output.Zbs;

Rwcmn  = Rbcmn(:,end);
Rwsmn  = Rbsmn(:,end);
Zwcmn  = Zbcmn(:,end);
Zwsmn  = Zbsmn(:,end);

% Compute (x,y) coordinates of the boundary surface

zeta   = zetaov2pi*(2*pi);

nth    = 2048;
dth    = 2*pi/nth;
theta  = dth:dth:2*pi; 

X      = zeros(1,nth);
Y      = zeros(1,nth);

for k=1:mn
 alpha  = double(im(k))*theta-double(in(k))*zeta;
 X = X + Rwcmn(k)*cos(alpha) + Rwsmn(k)*sin(alpha);
 Y = Y + Zwsmn(k)*sin(alpha) + Zwcmn(k)*cos(alpha);
end



% Plot Poincare section

if(newfig==1)
figure
end
hold on

scatter(X,Y,3,'filled', 'b')

axis equal
hold on
set(gca,'FontSize',12)
xlabel('R','FontSize',12)
ylabel('Z','FontSize',12)
%xlim([-1.1*rmax 1.1*rmax])
%ylim([-1.1*zmax 1.1*zmax])
