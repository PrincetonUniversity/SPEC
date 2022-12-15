function modify_Linitialize(filename, Linitialize)

sref  = ' Linitialize =';
snew1    = strcat(sref,{'      '}, num2str(Linitialize));
snew    = snew1{1};



%Open file 
fid   = fopen(filename,'rt');

tline = fgetl(fid);

count = 1;
lnum = 0;

A{1}  = tline;


% Read template file, copy lines in A, and identify reference lines

while ischar(tline)
    if(size(tline)>=size(sref))
      if(strcmp(tline(1:length(sref)),sref)==1)   
        lnum     = count;
      end
    end
  tline    = fgetl(fid);
  count    = count + 1;
  A{count} = tline;    
end

fclose(fid);


% Modify cell A at reference lines
A{lnum} = sprintf('%s',snew);

% Write cell A into new input file

fid = fopen(filename, 'w');

for i = 1:numel(A)
    if(A{i+1} == -1)
      fprintf(fid,'%s', A{i});
      break
    else
      fprintf(fid,'%s\n', A{i});
    end
end

fclose(fid);
end