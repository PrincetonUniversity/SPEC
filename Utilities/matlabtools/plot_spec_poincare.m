function plot_spec_poincare(data,nz0,nfp,trjstep,newfig)

% Produces Poincare plots of the field lines on different sections (within one field period)
%
% INPUT    
%   -data     : must be produced by calling read_spec(fname)
%   -nz0      : shows the nz0 toroidal plane or equidistant planes (nz0=-1)
%   -nfp      : number of field periods
%   -trjstep  : step to skip field-line trajectories when ploting (trjstep=0 means all trajectories are ploted)
%   -newfig   : opens(=1) or not(=0) a new figure, or overwrites (=2) last plot
%
% written by J.Loizu (2015)
% modified by A.Baillod (2019)
% modified by J.Loizu (2020)

try
 rpol    = data.input.physics.rpol; % get size of slab
catch
 rpol    = 1;
end

pdata    = pdata_from_data(data);   % construct stucture with poincare data information

nptraj   = size(pdata.R_lines,1);   % # of poincare trajectories (field lines)

nz       = size(pdata.R_lines,2);   % # of toroidal planes

flag2col = 'F';                     % flag for ploting field lines with alternating colour ('T') or not ('F')


disp(' ');
disp('Number of toroidal planes available: (one field period)');
nz
disp(' ');


rmax   = max(max(max(pdata.R_lines)));
rmin   = min(min(min(pdata.R_lines)));
zmax   = max(max(max(pdata.Z_lines)));
zmin   = min(min(min(pdata.Z_lines)));

switch pdata.Igeometry
    case 1
        xmin =  0;
        xmax =  2*pi*rpol;
        ymin = -0.1;
        ymax =  pdata.Rbc(1,end)+0.1;
    case 2
        xmin = -1.1*rmax;
        xmax =  1.1*rmax;
        ymin = -1.1*rmax;
        ymax =  1.1*rmax;
    case 3
        xmin =  0.9*rmin;
        xmax =  1.1*rmax;
        ymin =  1.1*zmin;
        ymax =  1.1*zmax;
end


nth    = 5096;  %ploting options for the boundary
bcol   = 'r';
bthick = 3;
if(pdata.Lfreebound==1)
bcol   = 'k';
bthick = 1;
end



if(flag2col=='T')
 pcol   = ['k' 'b'];
else
 pcol   = ['k' 'k'];
end

switch newfig
    case 0
        hold on
    case 1
        figure
        hold on
    case 2
        hold off
end

switch nz0

 case -1            

  npl = input('Number of planes: ');
  k   = 0;
  
  for j=1:nz/npl:nz  %for npl toroidal planes

   k = k+1;

   subplot(npl,1,k)
 
   switch pdata.Igeometry
       case 1
           R    = squeeze(pdata.R_lines(:,j,:));
           T    = rpol*mod(squeeze(pdata.th_lines(:,j,:)),2*pi);
           for i=1:1+trjstep:nptraj       %for each field line trajectory
               scatter(T(i,:),R(i,:),10,'.k')
               hold on
           end
       case 2  
           R    = squeeze(pdata.R_lines(:,j,:));
           T    = squeeze(pdata.th_lines(:,j,:));
           for i=1:1+trjstep:nptraj       %for each field line trajectory
               scatter(R(i,:).*cos(T(i,:)),R(i,:).*sin(T(i,:)),10,'.k')
               hold on;
           end
       case 3
           R = squeeze(pdata.R_lines(:,j,:));
           Z = squeeze(pdata.Z_lines(:,j,:));
           for i=1:1+trjstep:nptraj     %for each field line trajectory
            scatter(R(i,:),Z(i,:),10,'.k')
            hold on;
           end
       otherwise
           error('Unsupported geometry')
   end
   
   dth   = 2*pi/nth;
   theta = dth:dth:2*pi; 
   zeta  = (j-1)*(2*pi/nz)/nfp;

   switch pdata.Igeometry
     case 1
       Xb  = rpol*theta;
       Yb1 = 0;
       Yb2 = 0;
       for imn=1:pdata.mn     % get and plot the boundary
         alpha = double(pdata.im(imn))*theta-double(pdata.in(imn))*zeta;
         Yb1    = Yb1   + pdata.Rbc(imn,1)*cos(alpha) + pdata.Rbs(imn,1)*sin(alpha);
         Yb2    = Yb2   + pdata.Rbc(imn,end)*cos(alpha) + pdata.Rbs(imn,end)*sin(alpha);
       end
     case 2  
       Rb = 0;
       Zb = 0;
       for imn=1:pdata.mn     % get and plot the boundary
         alpha = double(pdata.im(imn))*theta-double(pdata.in(imn))*zeta;
         Rb = Rb + (pdata.Rbc(imn,end)*cos(alpha) + pdata.Rbs(imn,end)*sin(alpha)).*cos(theta);
         Zb = Zb + (pdata.Rbc(imn,end)*cos(alpha) + pdata.Rbs(imn,end)*sin(alpha)).*sin(theta);
       end
     case 3
       Rb = 0;
       Zb = 0;
       for imn=1:pdata.mn  % get and plot the boundary  % the pdata.in values go in steps of nfp
         alpha = double(pdata.im(imn))*theta-double(pdata.in(imn))*zeta;
         Rb    = Rb + pdata.Rbc(imn,end)*cos(alpha) + pdata.Rbs(imn,end)*sin(alpha);
         Zb    = Zb + pdata.Zbs(imn,end)*sin(alpha) + pdata.Zbc(imn,end)*cos(alpha); 
       end
     otherwise
       error('Unsupported geometry')
   end  

  if pdata.Igeometry ~= 1
   scatter(Rb,Zb,bthick,'*',bcol)
   hold on
   set(gca,'FontSize',12)
   axis equal
   xlabel('R','FontSize',12)
   ylabel('Z','FontSize',12)
  else
   scatter(Xb,Yb1,bthick,'*',bcol)
   scatter(Xb,Yb2,bthick,'*',bcol)
   hold on
   set(gca,'FontSize',12)
   xlabel('\theta r_{pol}','FontSize',12)
   ylabel('R','FontSize',12)
  end

  xlim([xmin xmax])
  ylim([ymin ymax])

 end

 otherwise  %if nz0>0
       
   switch pdata.Igeometry
       case 1
           R    = squeeze(pdata.R_lines(:,nz0,:));
           T    = rpol*mod(squeeze(pdata.th_lines(:,nz0,:)),2*pi);
           for i=1:1+trjstep:nptraj       %for each field line trajectory
               scatter(T(i,:),R(i,:),10,'.k')
               hold on
           end
       case 2  
           R    = squeeze(pdata.R_lines(:,nz0,:));
           T    = squeeze(pdata.th_lines(:,nz0,:));
           for i=1:1+trjstep:nptraj       %for each field line trajectory
               scatter(R(i,:).*cos(T(i,:)),R(i,:).*sin(T(i,:)),10,'.k')
               hold on;
           end
       case 3
           R = squeeze(pdata.R_lines(:,nz0,:));
           Z = squeeze(pdata.Z_lines(:,nz0,:));
           for i=1:1+trjstep:nptraj     %for each field line trajectory
            scatter(R(i,:),Z(i,:),10,'.k')
            hold on;
           end
       otherwise
           error('Unsupported geometry')
   end


  dth   = 2*pi/nth;
  theta = dth:dth:2*pi; 
  zeta  = (nz0-1.0)*(2.0*pi/nz)/double(nfp);
   
  
  switch pdata.Igeometry
     case 1
       Xb  = rpol*theta;
       Yb1 = 0;
       Yb2 = 0;
       for imn=1:pdata.mn     % get and plot the boundary
         alpha = double(pdata.im(imn))*theta-double(pdata.in(imn))*zeta;
         Yb1    = Yb1   + pdata.Rbc(imn,1)*cos(alpha) + pdata.Rbs(imn,1)*sin(alpha);
         Yb2    = Yb2   + pdata.Rbc(imn,end)*cos(alpha) + pdata.Rbs(imn,end)*sin(alpha);
       end
     case 2  
       Rb = 0;
       Zb = 0;
       for imn=1:pdata.mn     % get and plot the boundary
         alpha = double(pdata.im(imn))*theta-double(pdata.in(imn))*zeta;
         Rb = Rb + (pdata.Rbc(imn,end)*cos(alpha) + pdata.Rbs(imn,end)*sin(alpha)).*cos(theta);
         Zb = Zb + (pdata.Rbc(imn,end)*cos(alpha) + pdata.Rbs(imn,end)*sin(alpha)).*sin(theta);
       end
     case 3
       Rb = 0;
       Zb = 0;
       for imn=1:pdata.mn  % get and plot the boundary  % the pdata.in values go in steps of nfp
         alpha = double(pdata.im(imn))*theta-double(pdata.in(imn))*zeta;
         Rb    = Rb + pdata.Rbc(imn,end)*cos(alpha) + pdata.Rbs(imn,end)*sin(alpha);
         Zb    = Zb + pdata.Zbs(imn,end)*sin(alpha) + pdata.Zbc(imn,end)*cos(alpha); 
       end
     otherwise
       error('Unsupported geometry')
  end  

  if pdata.Igeometry ~= 1
   scatter(Rb,Zb,bthick,'*',bcol)
   hold on
   set(gca,'FontSize',12)
   axis equal
   xlabel('R','FontSize',12)
   ylabel('Z','FontSize',12)
  else
   scatter(Xb,Yb1,bthick,'*',bcol)
   scatter(Xb,Yb2,bthick,'*',bcol)
   hold on
   set(gca,'FontSize',12)
   xlabel('\theta r_{pol}','FontSize',12)
   ylabel('R','FontSize',12)
  end

  xlim([xmin xmax])
  ylim([ymin ymax])

end


disp(' ');
disp('--- end of program ---');
disp(' ');
