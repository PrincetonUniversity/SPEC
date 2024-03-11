function produce_spec_movie(inputroot,seqstart,seqstep,seqend,framerate,nfp,nz0,zetaov2pi,xrange,yrange,labtext)


%
% PRODUCE_SPEC_MOVIE( INPUTROOT, SEQSTART, SEQSTEP, SEQEND, FRAMERATE, NFP, NZ0, ZETAOV2PI, XRANGE, YRANGE, LABTEXT )
% ===================================================================================================================
%
% Produces Poincare movie from SPEC sequence of result files (Poincare data must exist)
%
% INPUT
% -----
%   -inputroot : spec input file name of the form 'somename_seq', and the '.end/.h5' files must already exist.
%   -seqstart  : sequence number to start the movie
%   -seqstep   : step to be taken as sequence files are read and showed
%   -seqend    : last sequence to be read and showed
%   -framerate : number of frames per second
%   -nfp       : Number of field periods 
%   -nz0       : toroidal plane number at which the Poincare plot is shown (e.g. nz0=1)
%   -zetao2pi  : toroidal angle over 2*pi at which the KAM surfaces are schown (zetao2pi=0 for nz0=1)
%   -xrange:    minimum and maximum values for xaxis, i.e. [xmin xmax]
%   -yrange:    minimum and maximum values for yaxis, i.e. [ymin ymax]
%   -labtext:   text label to add on each snapshot ('b' for beta, 'I' for current, 'none' for nothing)
%
%
%   written by J.Loizu (2017)
%   modified by J.Loizu (2019)

error('DEPRECATED: this needs to be reviewed')

%
writerObj           = VideoWriter('spec_movie.avi'); 

writerObj.FrameRate = framerate; %frames per second                 

open(writerObj); 
 
%
xmin = xrange(1);%8.5;%38.5;
xmax = xrange(2);%11.5;%41.5; 
ymin = yrange(1);%-1;%-1;
ymax = yrange(2);%1;%1; 
x1   = mean(xrange);%40;
y1   = 0.9*ymax;%0.9;%0.9;

mu0  = 1.25e-6; 
 
for it=seqstart:seqstep:seqend

 spec_hdf5   = strcat(inputroot,num2str(it),'.h5');

 data       = read_spec(spec_hdf5);

 if(labtext=='b')
 beta        = get_spec_beta(spec_hdf5);  
 figtext     = strcat('\beta_0 = ',num2str(round(beta*100/0.01)*0.01),' %');
 end
 
 if(labtext=='I')
 Itor        = get_spec_torcurr_kam_net(spec_hdf5);  
 figtext     = strcat('I_{\phi} = ', num2str(abs(round((Itor/mu0)/1000))),' kA');
 end
 
 plot_spec_poincare(data,nz0,nfp,0,1)
 
 plot_spec_kam(spec_hdf5,zetaov2pi,0)
 
 xlim([xmin xmax])
 ylim([ymin ymax])
  
 text(x1,y1,figtext);
 
 frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
 
 writeVideo(writerObj, frame);
 
 close
 
end

 close(writerObj); % Saves the movie.
