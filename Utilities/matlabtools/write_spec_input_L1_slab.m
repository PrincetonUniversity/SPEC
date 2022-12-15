function write_spec_input_L1_slab(template, inputname, Nvol, tfl, pfl, iota, oita, phit_edge, pre, mu, lrad, nptr)

% Writes spec input file from template with constraints corresponding to 
% Lconstraint=1 (tfl,iota,oita) for a slab geometry (where iota has size Nvol+1)
%
% INPUT
%   -template   : template input file name with .sp format 
%   -inputname  : new input file name with .sp format
%   -Nvol       : number of volumes
%   -tfl        : toroidal flux enclosed by each interface
%   -pfl	: poloidal flux enclosed by each interface
%   -iota       : rotational transform on the inner side of the interface
%   -oita       : rotational transform on the outer side of the interface
%   -phit_edge  : total toroidal flux
%   -pre        : pressure in each volume
%   -mu		: Lagrange multiplier in each volume
%   -lrad       : radial resolution in each volume
%   -nptr       : number of Poincare trajectories in each volume
%
%   written by J.Loizu (2019)

nlmod    = 10;

sref{1}  = ' pressure    =';
sref{2}  = ' tflux       =';
sref{3}  = ' iota        =';
sref{4}  = ' oita        =';
sref{5}  = ' Lrad        =';
sref{6}  = ' nPtrj       =';
sref{7}  = ' Nvol        =';
sref{8}  = ' pflux       =';
sref{9}  = ' mu          =';
sref{10} = ' phiedge     =';

snew1    = sref{1};
snew2    = sref{2};
snew3    = sref{3};%strcat(sref{3}, {'    '}, num2str(iota(1), 16), {'    '});
snew4    = sref{4};%strcat(sref{4}, {'    '}, num2str(oita(1), 16), {'    '});
snew5    = sref{5};
snew6    = sref{6};
snew7    = strcat(sref{7},{'   '},num2str(Nvol), {'   '});
snew8    = sref{8};
snew9    = sref{9};
snew10   = strcat(sref{10}, {'    '}, num2str(phit_edge), {'    '});

for i=1:Nvol
  snew1  = strcat(snew1, {'   '},num2str(pre(i) ,16), {'   '});
  snew2  = strcat(snew2, {'   '},num2str(tfl(i) ,16), {'   '});
  snew3  = strcat(snew3, {'   '},num2str(iota(i),16), {'   '});
  snew4  = strcat(snew4, {'   '},num2str(oita(i),16), {'   '});
  snew5  = strcat(snew5, {'   '},num2str(lrad(i),16), {'   '});
  snew6  = strcat(snew6, {'   '},num2str(nptr(i),16), {'   '});
  snew8  = strcat(snew8, {'   '},num2str(pfl(i) ,16), {'   '});
  snew9  = strcat(snew9, {'   '},num2str(mu(i)  ,16), {'   '});  
end

snew3  = strcat(snew3, {'   '},num2str(iota(Nvol+1),16), {'   '});
snew4  = strcat(snew4, {'   '},num2str(oita(Nvol+1),16), {'   '});


snew{1}    = snew1{1};
snew{2}    = snew2{1};
snew{3}    = snew3{1};
snew{4}    = snew4{1};
snew{5}    = snew5{1};
snew{6}    = snew6{1};
snew{7}    = snew7{1};
snew{8}    = snew8{1};
snew{9}    = snew9{1};
snew{10}   = snew10{1};



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
