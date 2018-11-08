function currents = run_spec_iter_fixtorcurr(execom,inputroot,Nvol,Lrad,nptr,startit,nit,Iref)


% Runs SPEC iteratively to find nearby equilibria with mu=0 in each volume and net toroidal surface-current Iref (units mu0*I).
% Target current is achieved by iterating on the enclosed poloidal flux.
%
% INPUT
%   - execom    : command for the execution of SPEC executable (e.g. './xspec' or 'mpirun -n 2 ./xspec')
%   - inputroot : spec input file name of the form 'somename_iter', and 'somename_iter0.sp.end/.h5' must already exist.
%   - Nvol      : number of volumes
%   - Lrad      : radial resolution in each volume (one value)
%   - nptr      : number of poincare trajectories traced in each volume
%   - startit   : starting iteration number
%   - nit       : number of iterations
%   - Iref      : targeted net-toroidal-current (array of size Nvol-1)
%
% OUTPUT
%   - currents  : array of net-toroidal-currents obtained on each interface at each interation step
%
% written by J.Loizu (2017)
% upgraded by J.Loizu (2018)


specexec = strcat(execom,{' '});

specexec = specexec{1};

currents = zeros(nit+1,Nvol-1);

pfac     = 1e-9;            % factor for the amplitude of the perturbation in pflux for Jacobian evaluation (reference: 1e-6)

ntheta   = 128;             % poloidal resolution for the loop integral evaluating the current

mu       = zeros(1,Nvol);   % as of now zero volume-currents are imposed


for i=1+startit:startit+nit

 display(' ');
 display(['START OF CURRENT ITERATION ' num2str(i)]);
 display(' ');

 %-------- Initial input and hdf5 files --------

 spec_input = strcat(inputroot,num2str(i-1),'.sp.end');
 spec_hdf5  = strcat(inputroot,num2str(i-1),'.sp.h5');

 tflux      = h5read(spec_hdf5,'/tflux');
 pflux      = h5read(spec_hdf5,'/pflux');
 pre        = h5read(spec_hdf5,'/pressure');
 IKAMtor    = get_spec_torcurr_kam_net(spec_hdf5,ntheta); 
 
 currents(i-startit,1:Nvol-1) = IKAMtor;
 
 %-------- Evaluation of J=dI/dX ---------

 tfl      = zeros(1,Nvol);
 
 pfl      = zeros(1,Nvol);

 dX       = zeros(1,Nvol);
 
 dIdX     = zeros(Nvol-1,Nvol-1);

 dX       = pflux*pfac;

 lrad     = Lrad*ones(1,Nvol);

 template = spec_input;


 % jacobian evaluation runs
   
 for lvol=2:Nvol
   
  display(' ');
  display(['JACOBIAN EVALUATION...' num2str(lvol-1)]);
  display(' ');
  
  tfl            = tflux; 

  pfl            = pflux; 
     
  pfl(lvol)      = pfl(lvol)  + dX(lvol);

  newinput       = strcat(inputroot,'jacobian',num2str(lvol),'.sp');

  write_spec_input_L0(template,newinput,Nvol,tfl,pfl,mu,pre,lrad,nptr);

  system(strcat([specexec newinput(1:end-3)]));
 
  spec_hdf5      = strcat(newinput(1:end-3),'.sp.h5');
  
  newIKAMtor     = get_spec_torcurr_kam_net(spec_hdf5,ntheta); 

  dIdX(:,lvol-1) = (newIKAMtor-IKAMtor)/dX(lvol);

 end
 
 
 %-------- Evaluation of dX steps from Newton method ---------

 Jinv                    = inv(dIdX);
 curr(1:Nvol-1)          = IKAMtor; 
 zfun(1:Nvol-1)          = IKAMtor-Iref; 
 
 dX(2:Nvol)              = -Jinv*transpose(zfun);
 
 display(' ');
 display(['CALCULATED STEP: ' num2str(transpose(dX)) ]);
 display(' ');
 
 %-------- Final run with appropriate X+dXstep ---------

 tfl            = tflux;
 pfl            = pflux;
 pfl(2:Nvol)    = pfl(2:Nvol) + dX(2:Nvol);

 newinput   = strcat(inputroot,num2str(i),'.sp');

 write_spec_input_L0(template,newinput,Nvol,tfl,pfl,mu,pre,lrad,nptr);

 system(strcat([specexec newinput(1:end-3)]));
 
 
 %-------- Update .DF file from last output --------

 system(strcat(['cp ' '.' inputroot num2str(i) '.sp.DF' ' .sp.DF']));

 
end

spec_hdf5                   = strcat(inputroot,num2str(startit+nit),'.sp.h5');
IKAMtor                     = get_spec_torcurr_kam_net(spec_hdf5,ntheta);
currents(nit+1,1:Nvol-1)    = IKAMtor;


