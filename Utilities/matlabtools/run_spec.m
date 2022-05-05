function run_spec(fname, comment)

% Go in ~/specruns/
!cd ~/specruns/

% Run SPEC
command = char(['./xspec ', fname]);
system(command)

% Move everything in there
h5fname = [fname, '.sp.h5'];

c = datetime('now');

matname = [fname, '.sp.mat'];


gdata = read_spec_grid(h5fname);
input = compose(string([   
            '\n Igeometry = ', num2str(gdata.Igeometry), ...
            '\n Lconstraint = ', num2str(gdata.Lconstraint), ...
            '\n Lfreebound = ', num2str(gdata.Lfreebound), ...
            '\n Nvol = ', num2str(gdata.Nvol), ...
            '\n Mrad = ', num2str(gdata.Mrad), ...
            '\n Nfp  = ', num2str(gdata.Nfp),  ...
            '\n Ntor = ', num2str(gdata.Ntor), ...
            '\n Mpol = ', num2str(gdata.Mpol), ...
            '\n']));
output = compose(string([  
            '\n ForceErr = ', num2str(gdata.ForceErr), ...
            '\n']));


out.date = c;
out.comment = comment;
out.input = fname;
out.h5name = h5fname;
out.input = input;
out.output = output;


% Write mat file
save(matname, '-struct', 'out');




end