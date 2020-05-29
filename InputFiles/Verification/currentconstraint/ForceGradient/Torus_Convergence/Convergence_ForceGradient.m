
init_fname = 'G3V08L3Fi.001.sp';

tmp = 1:1:15;
tmp = 1./2.^tmp;
DeltaR = tmp * 1E-1;

max_DeltaRel = zeros(1,length(DeltaR));

system(['../../../../../dspec ', init_fname])
hess = read_spec_hessian([init_fname, '.h5']);


% First do the runs
for ii = 1:length(DeltaR)
    
    dRZ = DeltaR(ii);

    fname_new = ['Run_', num2str(ii), '.sp'];
    change_spec_inputfile(init_fname, fname_new, 'dRZ', DeltaR(ii), 'Lcheck', 6, 'Lfindzero', 2)
    system(['../../../../../dspec ', fname_new])

   fname_out = ['Run_', num2str(ii) '_FD.txt'];
    system(['mv Lcheck6_output.FiniteDiff.txt ', fname_out])
    FG_FiniteDiff = importdata(fname_out);

    
    diff_abs_hess = abs(FG_FiniteDiff - hess);

    
    max_Delta_hess(ii) = max(max(diff_abs_hess));
end



figure
loglog(DeltaR, max_Delta_hess, '*')
ylabel('$log10(max(\mid\nabla F - \nabla F_{FD}\mid))$')
xlabel('$\log_{10}(\Delta R)$')
ax = gca;
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';


! rm Run* .Run* G3V08L3Fi.001.sp.* .G3V08L3Fi.001.sp*



