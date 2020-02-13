function plot_spec_poincare(data,nz0,nfp,trjstep,newfig)

% Produces Poincare plots of the field lines on different sections (within one field period)
%
% INPUT    
%   -data     : must be produced by calling read_spec_poincare(filename)
%   -nz0=-1   : shows a number of equidistant toroidal planes
%   -nz0=-2   : shows selected toroidal planes
%   -nz0>0    : shows the nz0 toroidal plane
%   -nfp      : is the number of field periods
%   -trjstep  : step to skip field-line trajectories when ploting (trjstep=0 means all trajectories are ploted)
%   -newfig   : opens(=1) or not(=0) a new figure. =2 to overwrite last
%   plot
%
%   originally written by J.Loizu (2015)


nptraj   = size(data.poincare.R,1);  % # of poincare trajectories (field lines)

nz       = size(data.poincare.R,2);  % # of toroidal planes

nppts    = size(data.poincare.R,3);  % # of iterations per trajectory

flag2col = 'F';                   % flag for ploting field lines with alternating colour ('T') or not ('F')


disp(' ');
disp('Number of toroidal planes available: (one field period)');
nz
disp(' ');


rmax   = max(max(max(data.poincare.R)));
rmin   = min(min(min(data.poincare.R)));
zmax   = max(max(max(data.poincare.Z)));
zmin   = min(min(min(data.poincare.Z)));

nth    = 5096;  %ploting options for the boundary
bcol   = 'r';
bthick = 3;
if(data.input.physics.Lfreebound==1)
bcol   = 'k';
bthick = 1;
end



if(flag2col=='T')
 pcol   = ['k' 'b'];
else
 pcol   = ['k' 'k'];
end

switch newfig			
	case 1			
    figure
    hold on;
	case 0
    hold on;
  case 2				% A.Baillod (06.2019)
    hold off;
end

switch nz0

 case -1            

  npl = input('Number of planes: ');
  k   = 0;
  
  for j=1:nz/npl:nz  %for npl toroidal planes

   k = k+1;

   subplot(npl,1,k)
 
   R = squeeze(data.poincare.R(:,j,:));

   Z = squeeze(data.poincare.Z(:,j,:));

   for i=1:nptraj     %for each field line trajectory
    scatter(R(i,:),Z(i,:),10,'.k')
   end

   Rb    = 0;    
   Zb    = 0;
   dth   = 2*pi/nth;
   theta = dth:dth:2*pi; 
   zeta  = (j-1)*(2*pi/nz)/nfp;

   for imn=1:data.output.mn  % get and plot the boundary  % the data.in values go in steps of nfp
    alpha = double(data.output.im(imn))*theta-double(data.output.in(imn))*zeta;
    Rb    = Rb + data.output.Rbc(imn,end)*cos(alpha) + data.output.Rbs(imn,end)*sin(alpha);
    Zb    = Zb + data.output.Zbs(imn,end)*sin(alpha) + data.output.Zbc(imn,end)*cos(alpha); 
   end

   scatter(Rb,Zb,bthick,'*',bcol)
   hold on;
    axis equal
    set(gca,'FontSize',12)
    xlim([0.9*rmin 1.1*rmax])
    ylim([1.1*zmin 1.1*zmax])
  end

 case -2            

  npl = input('Number of planes: ');
  ipl = ones(1,npl);

  for j=1:npl
   ipl(j) = input('plane number: ');
  end

  k   = 0;
  
  for j=1:npl %for npl toroidal planes

   k = k+1;

   subplot(npl,1,k)
 
   R = squeeze(data.poincare.R(:,ipl(j),:));

   Z = squeeze(data.poincare.Z(:,ipl(j),:));

   for i=1:nptraj     %for each field line trajectory
    scatter(R(i,:),Z(i,:),10,'.k')
   end

   Rb    = 0;    
   Zb    = 0;
   dth   = 2*pi/nth;
   theta = dth:dth:2*pi; 
   zeta  = (ipl(j)-1)*(2*pi/nz)/nfp;

   for imn=1:data.output.mn  % get and plot the boundary
    alpha = double(data.output.im(imn))*theta-double(data.output.in(imn))*zeta;
    Rb    = Rb + data.output.Rbc(imn,end)*cos(alpha) + data.output.Rbs(imn,end)*sin(alpha);
    Zb    = Zb + data.output.Zbs(imn,end)*sin(alpha) + data.output.Zbc(imn,end)*cos(alpha); 
   end

   scatter(Rb,Zb,bthick,'*',bcol)
   hold on;
    axis equal
    set(gca,'FontSize',12)
    xlim([0.9*rmin 1.1*rmax])
    ylim([1.1*zmin 1.1*zmax])
  end


 otherwise  
 
  R    = squeeze(data.poincare.R(:,nz0,:));

  Z    = squeeze(data.poincare.Z(:,nz0,:));

  for i=1:1+trjstep:nptraj       %for each field line trajectory 
   scatter(R(i,:),Z(i,:),10,'.',pcol(1+mod(i,2)))
   hold on
  end

  Rb    = 0;
  Zb    = 0;
  dth   = 2*pi/nth;
  theta = dth:dth:2*pi; 
  zeta  = (nz0-1.0)*(2.0*pi/nz)/double(nfp);			% modified by A.Baillod (06.2019)
  
  for imn=1:data.output.mn     % get and plot the boundary
   alpha = double(data.output.im(imn))*theta-double(data.output.in(imn))*zeta;
   Rb    = Rb + data.output.Rbc(imn,end)*cos(alpha) + data.output.Rbs(imn,end)*sin(alpha);
   Zb    = Zb + data.output.Zbs(imn,end)*sin(alpha) + data.output.Zbc(imn,end)*cos(alpha);
  end

  scatter(Rb,Zb,bthick,'*',bcol)
  
  axis equal
  hold on;
  set(gca,'FontSize',12)
  xlabel('R','FontSize',12)
  ylabel('Z','FontSize',12)
  xlim([0.9*rmin 1.1*rmax])
  ylim([1.1*zmin 1.1*zmax])
end


disp(' ');
disp('--- end of program ---');
disp(' ');
