function [x,y,z] = read_BIEST_vector(fname)

fid = fopen(fname, 'r');

dim = fread(fid, [2], 'uint64');
dim = dim(1);

data = fread(fid, [dim], 'double');

x = data(1:dim/3);
y = data(dim/3+1:2*dim/3);
z = data(2*dim/3+1:dim);

fclose(fid);

end

