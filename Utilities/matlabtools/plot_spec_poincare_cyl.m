function plot_spec_poincare_cyl(data,nz0,newfig)


% Produces Poincare plots of the field lines on different sections 
%   -data   must be produced by calling read_spec_poincare(filename)
%   -nz0=-1 shows a number of equidistant toroidal planes
%   -nz0=-2 shows selected toroidal planes
%   -nz0>0  shows the nz0 toroidal plane
%   -newfig opens(=1) or not(=0) a new figure
%   written by J.Loizu (2015)


nptraj = size(data.R_lines,1);  % # of poincare trajectories (field lines)

nz     = size(data.R_lines,2);  % # of toroidal planes

nppts  = size(data.R_lines,3);  % # of iterations per trajectory


disp(' ');
disp('Number of toroidal planes available:');
nz
disp(' ');

rmax   = max(max(max(data.R_lines)));

nth    = 5096;  %ploting options for the boundary
bcol   = 'r';
bthick = 10;

if(newfig==1)
figure
end
hold on

switch nz0

 case -1            

  npl = input(['Number of planes (divisor of ', num2str(nz), '): ']);
  k   = 0;
  
  for j=1:nz/npl:nz  %for npl toroidal planes

   k = k+1;

   subplot(npl,1,k)
 
   R = squeeze(data.R_lines(:,j,:));

   T = squeeze(data.th_lines(:,j,:));

   for i=1:nptraj     %for each field line trajectory
    scatter(R(i,:).*cos(T(i,:)),R(i,:).*sin(T(i,:)),10,'.k')
    axis equal
    hold on
    set(gca,'FontSize',12)
    xlabel('R','FontSize',12)
    ylabel('Z','FontSize',12)
    xlim([-1.1*rmax 1.1*rmax])
    ylim([-1.1*rmax 1.1*rmax])
   end

   xb    = 0;    
   yb    = 0;
   dth   = 2*pi/nth;
   theta = dth:dth:2*pi; 
   zeta  = (j-1)*(2*pi/nz);

   for imn=1:data.mn  % get and plot the boundary
    xb = xb + data.Rbc(imn,end)*cos(double(data.im(imn))*theta-double(data.in(imn))*zeta).*cos(theta);
    yb = yb + data.Rbc(imn,end)*cos(double(data.im(imn))*theta-double(data.in(imn))*zeta).*sin(theta);
   end

   scatter(xb,yb,bthick,'filled',bcol)

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
 
   R = squeeze(data.R_lines(:,ipl(j),:));

   T = squeeze(data.th_lines(:,ipl(j),:));

   for i=1:nptraj     %for each field line trajectory
    scatter(R(i,:).*cos(T(i,:)),R(i,:).*sin(T(i,:)),10,'.k')
    axis equal
    hold on
    set(gca,'FontSize',12)
    xlabel('R','FontSize',12)
    ylabel('Z','FontSize',12)
    xlim([-1.1*rmax 1.1*rmax])
    ylim([-1.1*rmax 1.1*rmax])
   end

   xb    = 0;    
   yb    = 0;
   dth   = 2*pi/nth;
   theta = dth:dth:2*pi; 
   zeta  = (ipl(j)-1)*(2*pi/nz);

   for imn=1:data.mn  % get and plot the boundary
    xb = xb + data.Rbc(imn,end)*cos(double(data.im(imn))*theta-double(data.in(imn))*zeta).*cos(theta);
    yb = yb + data.Rbc(imn,end)*cos(double(data.im(imn))*theta-double(data.in(imn))*zeta).*sin(theta);
   end

   scatter(xb,yb,bthick,'filled',bcol)

  end


 otherwise  
 
  R    = squeeze(data.R_lines(:,nz0,:));

  T    = squeeze(data.th_lines(:,nz0,:));

  for i=1:nptraj       %for each field line trajectory
   scatter(R(i,:).*cos(T(i,:)),R(i,:).*sin(T(i,:)),10,'.k')
   axis equal
   hold on
   set(gca,'FontSize',12)
   xlabel('R','FontSize',12)
   ylabel('Z','FontSize',12)
   xlim([-1.1*rmax 1.1*rmax])
   ylim([-1.1*rmax 1.1*rmax])
  end

  xb    = 0;
  yb    = 0;
  dth   = 2*pi/nth;
  theta = dth:dth:2*pi;
  zeta  = (nz0-1)*(2*pi/nz);

  for imn=1:data.mn     % get and plot the boundary
   xb = xb + data.Rbc(imn,end)*cos(double(data.im(imn))*theta-double(data.in(imn))*zeta).*cos(theta);
   yb = yb + data.Rbc(imn,end)*cos(double(data.im(imn))*theta-double(data.in(imn))*zeta).*sin(theta);
  end
  
  scatter(xb,yb,bthick,'filled',bcol)

end


disp(' ');
disp('--- end of program ---');
disp(' ');
