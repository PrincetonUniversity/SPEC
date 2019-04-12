function write_spec_input_vmec(template,inputname,specnl)

% Writes spec input file from template with data extracted from VMEC
%
% INPUT
%   -template   : template input file name with .sp format 
%   -inputname  : new input file name with .sp format
%   -specnl     : object contains things for replacement
%
%   written by J.Loizu (2016)
%   adapted by Zhisong Qu 12/2018

nlmod    = numel(specnl.sref);
format long
sref = specnl.sref;
for i = 1:nlmod
    sref{i}.string = sprintf(' %-12s=',sref{i}.name);
    sref{i}.newstring = sref{i}.string;
    for j = 1:numel(sref{i}.data)
        sref{i}.newstring = strcat(sref{i}.newstring, {'   '},num2str(sref{i}.data(j),16), {'   '});
    end
end
        
% Open template file for reading

fid   = fopen(template,'rt');

tline = fgetl(fid);

count = 1;

lnum  = zeros(1,nlmod);

A{1}  = tline;


% Read template file, copy lines in A, and identify reference lines
lastbs = 0;
rbcstart = 0;

while ischar(tline)
  for i=1:nlmod
    if(size(tline)>=size(sref{i}.string))
      if(strcmp(tline(1:length(sref{i}.string)),sref{i}.string)==1)   
        lnum(i)     = count;
      end      
    end
  end
  if length(tline) >= 1
      if tline == '/'
          lastbs = count;
      end
  end
  
  if length(tline) >= 3
      if strcmp(tline(1:3), 'Rbc') == 1
          lastrbc = count;
          count = count - 1;
      end
  end
  
  tline    = fgetl(fid);
  count    = count + 1;
  A{count} = tline;
end

fclose(fid);


% Modify cell A at reference lines

for i=1:nlmod
  A{lnum(i)} = sprintf('%s',sref{i}.newstring{1});
end

% Write cell A into new input file

fid = fopen(inputname, 'w');

for i = 1:lastbs
    
    if i == lastrbc 
        for j = 1:specnl.interfaces.mn
            for k = 1:numel(specnl.sref2d)
              fprintf(fid,'%s(%3d,%3d)=%23.16e, ', specnl.sref2d{k}.name, specnl.interfaces.in(j), ...
                  specnl.interfaces.im(j), specnl.sref2d{k}.data(j));
            end
            fprintf(fid, '\n');
        end
    end
    
    if(A{i+1} == -1)
      fprintf(fid,'%s', A{i});
      break
    else
      fprintf(fid,'%s\n', A{i});
    end
end

% Write initial guess of the interfaces
for i = 1:specnl.interfaces.mn
    fprintf(fid, '%3d %3d ', specnl.interfaces.im(i), specnl.interfaces.in(i));
    for j = 2:specnl.interfaces.nvol+1
        fprintf(fid, '%23.16e %23.16e %23.16e %23.16e ', specnl.interfaces.irbc(i,j), ...
            specnl.interfaces.izbs(i,j), -specnl.interfaces.irbs(i,j), ...
            specnl.interfaces.izbc(i,j));
    end
    fprintf(fid, '\n');
end

fclose(fid);


