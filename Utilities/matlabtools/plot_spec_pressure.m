function plot_spec_pressure(fname)

% Plots stepped-pressure profile (without pscale) versus normalized toroidal flux used in SPEC
%   -fname   : filename in HDF5 format 
%   written by J.Loizu (2018)


data = read_spec_grid(fname);

pvol = data.pressure;

tfl  = data.tflux;

Nvol = data.Nvol;

p0   = zeros(1,10);

figure
hold on

p0(1:end) = pvol(1);
tmin      = 0;
tmax      = tfl(1);
tarr      = linspace(tmin,tmax,10);
plot(tarr,p0,'b')
x         = [tfl(1),tfl(1)];
y         = [pvol(2) pvol(1)];
plot(x,y,'b')

for i=2:Nvol-1

 p0(1:end) = pvol(i);
 tmin      = tfl(i-1);
 tmax      = tfl(i);
 tarr      = linspace(tmin,tmax,10);
 plot(tarr,p0,'b')
 x         = [tfl(i),tfl(i)];
 y         = [pvol(i+1) pvol(i)];
 plot(x,y,'b')

end

p0(1:end) = pvol(Nvol);
tmin      = tfl(Nvol-1);
tmax      = tfl(Nvol);
tarr      = linspace(tmin,tmax,10);
plot(tarr,p0,'b')
x         = [tfl(Nvol),tfl(Nvol)];
y         = [0 pvol(Nvol)];
plot(x,y,'b')

ylabel('p')
xlabel('\Psi / \Psi_{edge}')



