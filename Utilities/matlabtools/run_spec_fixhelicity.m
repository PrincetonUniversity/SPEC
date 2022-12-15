function helicities = run_spec_fixhelicity(flagloc,inputroot,Nvol,Lrad,nptr,startseq,nit,Href)


% Runs SPEC iteratively in order to find nearby equilibria with given toroidal and poloidal fluxes, AND helicity in each volume.
% Toroidal and poloidal fluxes are enforced via Lconstraint=0. Target helicity is achieved by iterating on the beltrami parameter mu.
%   -flagloc   : flag to specify whether the SPEC is run locally (1) or on the cluster (0)
%   -inputroot : spec input file name of the form 'somename_seq', and 'somename_seq0.sp.end/.h5' must already exist.
%   -Nvol      : number of volumes
%   -Lrad      : radial resolution in each volume
%   -nptr      : number of poincare trajectories traced in each volume
%   -startseq  : starting sequence number
%   -nit       : number of iterations
%   -Href      : targeted helicity (array of size Nvol)
%
% Note: so far written for slab geometry
%
%   written by J.Loizu (2018)


pfac = 1e-6; % factor for the amplitude of the perturbation in pflux for Jacobian evaluation (ref: 1e-6)


if(flagloc==1)
 specexec = './xspec ';
else
 specexec = 'mpirun ./xspec_clus ';
end

helicities  = zeros(nit+1,Nvol);


for i=1+startseq:startseq+nit

 display(' ');
 display(['START OF HELICITY ITERATION ' num2str(i)]);
 display(' ');

 %-------- Initial input and hdf5 files --------

 spec_input = strcat(inputroot,num2str(i-1),'.sp.end');
 spec_hdf5  = strcat(inputroot,num2str(i-1),'.sp.h5');

 tflux      = h5read(spec_hdf5,'/tflux');
 pflux      = h5read(spec_hdf5,'/pflux');
 mu         = h5read(spec_hdf5,'/mu');
 pre        = h5read(spec_hdf5,'/pressure');
 hel        = h5read(spec_hdf5,'/helicity');
 
 helicities(i-startseq,1:Nvol) = hel;
 
 %-------- Evaluation of J=dH/dX ---------

 tfl      = zeros(1,Nvol);
 
 pfl      = zeros(1,Nvol);

 mul      = zeros(1,Nvol);

 dX       = zeros(1,Nvol);
 
 dHdX     = zeros(Nvol,Nvol);

 dX       = max(abs(mu))*ones(Nvol,1)*pfac;     % max(mu) to avoid zeroes when mu(lvol)=0

 lrad     = Lrad*ones(1,Nvol);

 template = spec_input;


 % jacobian evaluation runs
   
 for lvol=1:Nvol
   
  display(' ');
  display(['JACOBIAN EVALUATION...' num2str(lvol)]);
  display(' ');
  
  tfl            = tflux; 

  pfl            = pflux; 

  mul            = mu;
     
  mul(lvol)      = mul(lvol)  + dX(lvol);

  newinput       = strcat(inputroot,'jacobian',num2str(lvol),'.sp');

  write_spec_input_L0(template,newinput,Nvol,tfl,pfl,mul,pre,lrad,nptr);

  system(strcat([specexec newinput(1:end-3)]));
 
  spec_hdf5      = strcat(newinput(1:end-3),'.sp.h5');
  
  newhel         = h5read(spec_hdf5,'/helicity');

  dHdX(:,lvol)   = (newhel-hel)/dX(lvol);

 end
 
 
 %-------- Evaluation of dX steps from Newton method ---------

 Jinv                    = inv(dHdX);
 zfun(1:Nvol)            = hel-Href;
 
 dX(1:Nvol)              = -Jinv*transpose(zfun);
 
 display(' ');
 display(['CALCULATED STEP: ' num2str(transpose(dX)) ]);
 display(' ');
 
 %-------- Final run with appropriate X+dXstep ---------

 tfl            = tflux;
 pfl            = pflux;
 mul            = mu;
 mul(1:Nvol)    = mul(1:Nvol) + dX(1:Nvol);

 newinput   = strcat(inputroot,num2str(i),'.sp');

 write_spec_input_L0(template,newinput,Nvol,tfl,pfl,mul,pre,lrad,nptr);

 system(strcat([specexec newinput(1:end-3)]));
 
 
 %-------- Update .GF file from last output --------

% system(strcat(['cp ' '.' inputroot num2str(i) '.sp.DF' ' .DF']));
% system(strcat(['cp ' '.' inputroot num2str(i) '.sp.DF' ' ' inputroot num2str(i+1) '.sp.DF']));
 
end

spec_hdf5                   = strcat(inputroot,num2str(startseq+nit),'.sp.h5');
hel                         = h5read(spec_hdf5,'/helicity');
helicities(nit+1,1:Nvol)    = hel;


