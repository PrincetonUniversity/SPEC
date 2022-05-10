function change_spec_inputfile(fname, fname_new, varargin)

% CHANGE_SPEC_INPUTFILE
% ---------------------
%
% Generate a new input file for SPEC from a template, where any number of
% parameter can be changed
%
% INPUT
% -----
%   fname:      Name of the template. This file will not be modified by the
%               execution of the script
%   fname_new:  Name of the new file. This file will be generated from the
%               template with the required modifications
%   modifs:     Any pair of field - value present in the input file. For
%               example 'Nvol', 8, 'Lconstraint', 1 will change the number 
%               of volumes to 8 and the constraint to 1.
%
% OUTPUT
% ------
%   A new SPEC input file in the directory of execution
%
% Written by A. Baillod (2020)

Nfield = length(varargin);
if mod(Nfield,2)~=0
    error('Incorrect input list')
end

Nfield = Nfield /2;

sref = cell(1, Nfield);
snew = cell(1, Nfield);
for jj=1:Nfield
    
    field = varargin(2*jj -1);
    value = varargin(2*jj);
    
    field = field{1}; value = value{1};
    
    N = length(value);

    sref{jj}  = [' ',field];

    snew{jj}  = sref{jj};
    nsnew = length(snew{jj});
    diff = 14 - nsnew;
    for ii=1:diff
       snew{jj} = [snew{jj}, ' ']; 
    end

    snew{jj} = [snew{jj}, '= '];

    for ii=1:N
       if isa(value(ii), 'double')
        snew{jj} = strcat(snew{jj}, {'   '},num2str(value(ii) ,16), {'   '});  
       elseif isa(value(ii), 'char')
        snew{jj} = strcat(snew{jj}, {'   '}, value(ii)            , {'   '}); 
       else
        error('Unsupported data type')
       end
           
    end

    snew{jj} = snew{jj}{1};
end
    
% Open template file for reading

fid   = fopen(fname,'rt');

tline = fgetl(fid);

count = 1;

lnum  = zeros(1,Nfield);

A{1}  = tline;

% Read template file, copy lines in A, and identify reference lines

while ischar(tline)
  for i=1:Nfield
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

for i=1:Nfield
  A{lnum(i)} = sprintf('%s',snew{i});
end


% Write cell A into new input file

fid = fopen(fname_new, 'w');

for i = 1:numel(A)
    if(A{i+1} == -1)
      fprintf(fid,'%s', A{i});
      break
    else
      fprintf(fid,'%s\n', A{i});
    end
end

fclose(fid);
