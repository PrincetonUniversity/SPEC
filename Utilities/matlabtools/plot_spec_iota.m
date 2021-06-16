function out = plot_spec_iota(data,iorq,xaxis,newfig)

%
% PLOT_SPEC_IOTA( DATA, IORQ, XAXIS, NEWFIG )
% ===========================================
%
% Produces rotational transform plot from field-line-tracing data
%
% INPUT
% -----
%   -data      : produced by calling read_spec(fname)
%   -iorq      : plot iota('i') or safety factor ('q')
%   -xaxis='s' : plots s-coordinnate as the x-axis
%   -xaxis='R' : plots R-coordinnate as the x-axis
%   -xaxis='f' : plots toroidal flux as the x-axis
%   -xaxis='r' : plots sqrt(toroidal flux) as the x-axis
%   -newfig    : opens(=1) or not(=0) a new figure, or overwrites(=2) current plot
%
% OUTPUT
% ------
%   -out       : structure with arrays of iota and the chosen xaxis array
%
% written by J.Loizu (2015)
% modified by J.Loizu (06.2017)
% debugged by J.Loizu (02.2018)
% modified by A.Baillod (01.2019)
% modified by J.Loizu (01.2020)

if(newfig==1)
    figure
    hold on;
elseif newfig==0
    hold on;
elseif newfig==2			% A.Baillod (06.2019)
    hold off;
end

if(iorq=='i')
F = data.transform.fiota(1:end,2);
Flabel='\iota';
elseif(iorq=='q')
F = 1./data.transform.fiota(1:end,2);
Flabel='q';
end

nsucctrj = length(data.poincare.R(:,1,1)); % number of successfully followed trajectories

switch xaxis
 case 's'
  plot(data.transform.fiota(1:nsucctrj,1),F(1:nsucctrj),'*','MarkerSize',8,'LineWidth',2)
  ylabel(Flabel)
  out    = cell(2);
  out{1} = data.transform.fiota(:,1); 
  out{2} = F;

 case 'R'
  plot(transpose(data.poincare.R(:,1,1)),F(1:nsucctrj),'*','MarkerSize',8,'LineWidth',2)
  ylabel(Flabel)
  out = cell(2);
  out{1} = data.poincare.R(:,1,1);            
  out{2} = F(1:nsucctrj);  

 case 'f'
  sval    = data.transform.fiota(1:nsucctrj,1);
  nvol    = data.input.physics.Nvol;
  mvol    = data.output.Mvol;
  nptrj   = zeros(1,mvol);
  count   = 1;
  
  ind = find(sval==-2); %Remove wrongly written data
  sval(ind) = [];
  F(ind) = [];
  nsucctrj = nsucctrj - length(ind);
  
  for is=1:length(sval)-1
   if(sval(is)>0 && sval(is+1)<0) %Reached the end of a volume
    if(count==1) 
    nptrj(count) = is;
    else
    nptrj(count) = is-sum(nptrj);
    end
    count        = count+1;
   end
  end
  nptrj(mvol)    = length(sval)-sum(nptrj);
  
  ns       = 32;
  nt       = 32;
  cumflux  = 0;
  kstart   = 1;
  psitor   = zeros(1,length(sval));


  for lvol=1:mvol
      
      if lvol==nvol % Required if last value had a -2 (often in free bound)
          phiedge = cumflux + get_spec_torflux(data, nvol, 0, -1, 1, ns, nt);
      end
      
      for k=kstart:kstart-1+nptrj(lvol)
          psitor(k) = cumflux + get_spec_torflux(data,lvol,0,-1,sval(k),ns,nt);
      end
      cumflux = psitor(k);
      kstart  = kstart+nptrj(lvol);
  end
    
  plot(psitor/phiedge,F(1:nsucctrj),'*','MarkerSize',8,'LineWidth',2)
  ylabel(Flabel)
  xlabel('\Psi / \Psi_{edge}')
  
  out    = cell(2);
  out{1} = psitor/psitor(end);
  out{2} = F(1:nsucctrj);
  
case 'r'
  sval    = data.transform.fiota(1:nsucctrj,1);
  nvol    = data.input.physics.Nvol;
  mvol    = data.output.Mvol;
  nptrj   = zeros(1,nvol);
  count   = 1;
  
  ind = find(sval==-2); %Remove wrongly written data
  sval(ind) = [];
  F(ind) = [];
  nsucctrj = nsucctrj - length(ind);
  
  for is=1:length(sval)-1
   if(sval(is)>0 && sval(is+1)<0)
    if(count==1) 
    nptrj(count) = is;
    else
    nptrj(count) = is-sum(nptrj);
    end
    count        = count+1;
   end
  end
  nptrj(mvol)    = length(sval)-sum(nptrj);

  ns       = 32;
  nt       = 32;
  cumflux  = 0;
  kstart   = 1;
  psitor   = zeros(1,length(sval));

  for lvol=1:mvol
      
      if lvol==nvol % Required if last value had a -2 (often in free bound)
          phiedge = cumflux + get_spec_torflux(data, nvol, 0, -1, 1, ns, nt);
      end
      
      for k=kstart:kstart-1+nptrj(lvol)
          psitor(k) = cumflux + get_spec_torflux(data,lvol,0,-1,sval(k),ns,nt);
      end
      cumflux = psitor(k);
      kstart  = kstart+nptrj(lvol);
  end
   
%  phiedge = data.input.physics.phiedge;

  plot(sqrt(psitor/phiedge),F(1:nsucctrj),'*','MarkerSize',8,'LineWidth',2)
  ylabel(Flabel)
  xlabel('(\Psi / \Psi_{edge})^{1/2}')
  
  out    = cell(2);
  out{1} = psitor/psitor(end);
  out{2} = F(1:nsucctrj);
  
end
