function plot_spec_kam_slab(filename,zetaov2pi,newfig)

% Produces a "Poincare plot" of the ideal-interfaces in slab geometry
%
% INPUT
%   -filename   : path to the hdf5 output file (e.g. 'testcase.sp.h5')
%   -zetaov2pi  : shows the toroidal plane at zeta=2*pi*(zetaov2pi) 
%   -newfig     : opens(=1) or not(=0) a new figure
%
%   written by J.Loizu (2017)


Nvol   = h5read(filename,'/Nvol');
tflux  = h5read(filename,'/tflux');
mn     = h5read(filename,'/mn');
im     = h5read(filename,'/im');
in     = h5read(filename,'/in');
Rbcmn  = h5read(filename,'/Rbc');
Rbsmn  = h5read(filename,'/Rbs');
try
 rpol   = h5read(filename,'/rpol');
catch
 rpol   = 1;
end

% Compute (x,y) coordinates of each interface

zeta   = zetaov2pi*(2*pi);
nth    = 2048;
dth    = 2*pi/nth;
theta  = dth:dth:2*pi; 

X      = zeros(Nvol,nth);
Y      = zeros(Nvol,nth);

for i=1:Nvol
 X(i,:) = theta*rpol;
 for k=1:mn
 alpha  = double(im(k))*theta-double(in(k))*zeta;
 Y(i,:) = Y(i,:) + Rbcmn(k,i+1)*cos(alpha) + Rbsmn(k,i+1)*sin(alpha);
 end
end


% Plot "Poincare section"

if(newfig==1)
figure
end
hold on

for i=1:Nvol
 scatter(X(i,:),Y(i,:),3,'filled','r')
end

set(gca,'FontSize',12)
xlabel('\theta r_{pol}','FontSize',12)
ylabel('R','FontSize',12)
xlim([min(min(X)) max(max(X))])
ylim([min(min(Y)) max(max(Y))])

