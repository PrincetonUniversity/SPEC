function plot_spec_kam_slab(filename,zetaov2pi,newfig)

% Produces a "Poincare plot" of the ideal-interfaces in slab geometry
%
% INPUT
%   -filename  : path to the hdf5 output file (e.g. 'testcase.sp.h5')
%   -nz0       : the toroidal plane to be shown
%   -newfig    : opens(=1) or not(=0) a new figure.
%
%   written by J.Loizu (2017)


Nvol   = h5read(filename,'/Nvol');
tflux  = h5read(filename,'/tflux');
mn     = h5read(filename,'/mn');
im     = h5read(filename,'/im');
in     = h5read(filename,'/in');
Rbcmn  = h5read(filename,'/Rbc');
Rbsmn  = h5read(filename,'/Rbs');


% Compute (x,y) coordinates of each interface

zeta   = zetaov2pi*(2*pi);
nth    = 2048;
dth    = 2*pi/nth;
theta  = dth:dth:2*pi; 

X      = zeros(Nvol,nth);
Y      = zeros(Nvol,nth);

for i=1:Nvol
 X(i,:) = theta;
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


