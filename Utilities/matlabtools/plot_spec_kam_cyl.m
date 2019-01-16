function plot_spec_kam_cyl(filename,nz0,newfig)


% Produces a "Poincare plot" of the KAM surfaces in cylindrical geometry
%   -filename   SPEC output hdf5 file
%   -nz0        shows the nz0 toroidal plane
%   -newfig     opens(=1) or not(=0) a new figure
%   written by J.Loizu (2015)

period = 1;

Nvol   = h5read(filename,'/Nvol');
tflux  = h5read(filename,'/tflux');
mn     = h5read(filename,'/mn');
im     = h5read(filename,'/im');
in     = h5read(filename,'/in');
Rmn    = h5read(filename,'/Rbc');

rmax   = Rmn(1,end);


% Compute (x,y) coordinates of each KAM surface

zeta   = 0*pi;

nth    = 2048;
dth    = 2*pi/nth;
theta  = dth:dth:2*pi; 

X      = zeros(Nvol,nth);
Y      = zeros(Nvol,nth);

for i=1:Nvol
 for k=1:mn
 X(i,:) = X(i,:) + Rmn(k,i+1)*cos(double(im(k))*theta-double(in(k))*zeta).*cos(theta);
 Y(i,:) = Y(i,:) + Rmn(k,i+1)*cos(double(im(k))*theta-double(in(k))*zeta).*sin(theta);
 end
end


% Plot Poincare section

if(newfig==1)
figure
end
hold on

for i=1:double(Nvol)
 if(double(i/period-floor(i/period))==0)
 scatter(X(i,:),Y(i,:),3,'filled','r')
 end
end
axis equal
hold on
set(gca,'FontSize',12)
xlabel('R','FontSize',12)
ylabel('Z','FontSize',12)
xlim([-1.1*rmax 1.1*rmax])
ylim([-1.1*rmax 1.1*rmax])
