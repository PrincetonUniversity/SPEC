function plot_spec_iotakam_git(filename,freebound,iorq,xaxis,newfig)


% Plots rotational transform on each side of each ideal interface
%
% INPUT
%   -filename   : path to the hdf5 output file (e.g. 'testcase.sp.h5')
%   -freebound  : flag for whether it is a fixed (0) or free (1) boundary calculation
%   -iorq       : plot iota('i') or safety factor ('q')
%   -xaxis='R'  : plots R-position of interfaces (at phi=0) as the x-axis
%   -xaxis='f'  : plots toroidal flux as the x-axis
%   -xaxis='r'  : plots sqrt(toroidal flux) as the x-axis
%   -newfig     : opens(=1) or not(=0) a new figure
%
% written by J.Loizu (2017) 



Nvol   = h5read(filename,'/Nvol');
tflux  = h5read(filename,'/tflux');
iota   = h5read(filename,'/iota');
oita   = h5read(filename,'/oita');
Rmn    = h5read(filename,'/Rbc');
im     = h5read(filename,'/im');
in     = h5read(filename,'/in');
mn     = h5read(filename,'/mn');

R0     = zeros(1,Nvol);

for l=1:Nvol
 R0(l) = sum(Rmn(:,l+1));
end

if(newfig==1)
figure
end
hold on

if(iorq=='i')
F = abs(iota(2:end));
G = abs(oita(2:end));
Flabel='\iota';
elseif(iorq=='q')
F = 1./abs(iota(2:end));
G = 1./abs(oita(2:end));
Flabel='q';
end


switch xaxis
 case 'R'
  plot(R0,F,'r+','MarkerSize',6,'LineWidth',2)
  plot(R0,G,'m+','MarkerSize',6,'LineWidth',2)
  xlabel('R')
  ylabel(Flabel)
 case 'f'
  plot(tflux(1:end-freebound),F,'r+','MarkerSize',6,'LineWidth',2)
  plot(tflux(1:end-freebound),G,'m+','MarkerSize',6,'LineWidth',2)
  xlabel('\Psi / \Psi_{edge}')
  ylabel(Flabel)
 case 'r'
  plot(sqrt(tflux(1:end-freebound)),F,'r+','MarkerSize',6,'LineWidth',2)
  plot(sqrt(tflux(1:end-freebound)),G,'m+','MarkerSize',6,'LineWidth',2)
  xlabel('(\Psi / \Psi_{edge})^{1/2}')
  ylabel(Flabel)  
end
