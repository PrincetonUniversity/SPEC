function write_spec_input_L3(template, inputname, Nvol, Lfreebound, tfl, phit_edge, ...
                             pflux, Ivol, Isurf, curtor, pressure, lrad, nptr, Linitialize)

% 
% WRITE_SPEC_INPUT_L3( TEMPLATE, INPUTNAME, NVOL, LFREEBOUND, TFL, PHIT_EDGE, ... )
% =================================================================================
%
% Writes spec input file from template with constraints corresponding to 
% Lconstraint=1 (tfl,iota,oita)
%
% INPUT
% -----
%   -template   : template input file name with .sp format 
%   -inputname  : new input file name with .sp format
%   -Nvol       : number of volume
%   -Lfreebound : Flag for free boundary
%   -tfl        : toroidal flux enclosed by each interface
%   -phit_edge  : total toroidal flux
%   -pflux      : Poloidal flux (used as guess)
%   -Ivol       : volume current in each volume
%   -Isurf      : surface current at each interface
%   -curtor     : Total plasma current. only relevant if Lfreebound=1
%   -pressure   : pressure in each volume
%   -lrad       : radial resolution in each volume
%   -nptr       : number of Poincare trajectories in each volume
%   -Linitialize: How to initialize interfaces
%
%   written by A.Baillod (2019)

Mvol = Nvol+Lfreebound;

nlmod    = 11;

sref{1}  = ' pressure    =';
sref{2}  = ' tflux       =';
sref{3}  = ' Lrad        =';
sref{4}  = ' nPtrj       =';
sref{5}  = ' Nvol        =';
sref{6}  = ' phiedge     =';
sref{7}  = ' Ivolume     =';
sref{8}  = ' Isurf       =';
sref{9}  = ' curtor      =';
sref{10}  = ' pflux       =';
sref{11}  = ' Linitialize =';

snew1    = sref{1};
snew2    = sref{2};
snew3    = sref{3};
snew4    = sref{4};
snew5    = strcat(sref{5}, {'    '}, num2str(Nvol)     , {'    '});
snew6    = strcat(sref{6}, {'    '}, num2str(phit_edge), {'    '});
snew7    = sref{7};
snew8    = sref{8};
snew9    = strcat(sref{9}, {'    '}, num2str(curtor), {'    '});
snew10   = sref{10};
snew11    = strcat(sref{11}, {'    '}, num2str(Linitialize), {'    '});

for i=1:Mvol
  snew1  = strcat(snew1, {'   '},num2str(pressure(i) ,16), {'   '});
  snew2  = strcat(snew2, {'   '},num2str(tfl(i)      ,16), {'   '});
  snew3  = strcat(snew3, {'   '},num2str(lrad(i)     ,16), {'   '});
  snew4  = strcat(snew4, {'   '},num2str(nptr(i)     ,16), {'   '});
  snew7  = strcat(snew7, {'   '},num2str(Ivol(i)     ,16), {'   '});
  snew8  = strcat(snew8, {'   '},num2str(Isurf(i)    ,16), {'   '}); 
  snew10 = strcat(snew10, {'   '},num2str(pflux(i)    ,16), {'   '}); 
end

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
snew{11}   = snew11{1};


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
