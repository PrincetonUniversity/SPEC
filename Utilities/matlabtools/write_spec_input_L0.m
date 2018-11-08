function write_spec_input_L0(template,inputname,Nvol,tfl,pfl,mu,pre,lrad,nptr)

% Writes spec input file from template with constraints corresponding to Lconstraint=0 (tfl,pfl,mu)
%
% INPUT
%   -template   : template input file name with .sp format 
%   -inputname  : new input file name with .sp format
%   -Nvol       : number of volumes
%   -tfl        : toroidal flux enclosed by each interface
%   -pfl        : poloidal flux enclosed by each interface
%   -mu         : beltrami parameter in each volume
%   -pre        : pressure in each volume
%   -lrad       : radial resolution in each volume
%   -nptr       : number of Poincare trajectories in each volume
%
%   written by J.Loizu (2016)

nlmod    = 6;

sref{1}  = ' pressure    =';
sref{2}  = ' tflux       =';
sref{3}  = ' pflux       =';
sref{4}  = ' mu          =';
sref{5}  = ' Lrad        =';
sref{6}  = ' nPtrj       =';

snew1    = sref{1};
snew2    = sref{2};
snew3    = sref{3};
snew4    = sref{4};
snew5    = sref{5};
snew6    = sref{6};

for i=1:Nvol
  snew1  = strcat(snew1,{'   '},num2str(pre(i),16), {'   '});
  snew2  = strcat(snew2,{'   '},num2str(tfl(i),16), {'   '});
  snew3  = strcat(snew3,{'   '},num2str(pfl(i),16), {'   '});
  snew4  = strcat(snew4,{'   '},num2str(mu(i),16),  {'   '});
  snew5  = strcat(snew5,{'   '},num2str(lrad(i),16),{'   '});
  snew6  = strcat(snew6,{'   '},num2str(nptr(i),16),{'   '});
end

snew{1}    = snew1{1};
snew{2}    = snew2{1};
snew{3}    = snew3{1};
snew{4}    = snew4{1};
snew{5}    = snew5{1};
snew{6}    = snew6{1};

% Open template file for reading

fid   = fopen(template,'rt');

tline = fgetl(fid);

count = 1;

lnum  = zeros(1,nlmod);

A{1}  = tline;


% Read template file, copy lines in A, and identify reference lines

while ischar(tline)
  for i=1:nlmod
    if(size(tline)>=size(sref{i}))
      if(strcmp(tline(1:length(sref{i})),sref{i})==1)   
        lnum(i)     = count;
      end
    end
  end
  tline    = fgetl(fid);
  count    = count + 1;
  A{count} = tline;    
end

fclose(fid);


% Modify cell A at reference lines

for i=1:nlmod
  A{lnum(i)} = sprintf('%s',snew{i});
end

% Write cell A into new input file

fid = fopen(inputname, 'w');

for i = 1:numel(A)
    if(A{i+1} == -1)
      fprintf(fid,'%s', A{i});
      break
    else
      fprintf(fid,'%s\n', A{i});
    end
end

fclose(fid);


