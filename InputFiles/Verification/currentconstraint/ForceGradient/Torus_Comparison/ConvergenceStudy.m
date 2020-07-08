init_fname = 'G3V08L3Fr.001.sp';


max_DeltaRel = [];

tmp = 1:15;
tmp = 1./2.^tmp;
DeltaR = tmp * 1E-2;

% First do the runs
for ii = 1:length(DeltaR)
    
    dRZ = DeltaR(ii);

    fname_new = ['Run_', num2str(ii), '.sp'];
    change_spec_inputfile(init_fname, fname_new, 'dRZ', DeltaR(ii), 'Lcheck', 6, 'Lfindzero', 2)
    system(['../../../../../dspec ', fname_new])

    fname_out = ['Run_', num2str(ii) '.Lcheck6_output.FiniteDiff.txt'];
    FG_FiniteDiff = importdata(fname_out);

    FG_analytical = importdata('Run_1.Lcheck6_output.txt');

    diff_abs = abs(FG_analytical - FG_FiniteDiff);

    max_DeltaRel(ii) = max(max(diff_abs)) / max(max(abs(FG_analytical)));
end


figure('Position', [200,200,900,700])
grid on;
ax = loglog(DeltaR, max_DeltaRel, '*', 'LineWidth', 2.3, 'MarkerSize', 12);
xlabel('$\Delta R$')
ylabel('$\max\mid\nabla F - \nabla F_{FD}\mid / \max\mid\nabla F\mid$')
ax = gca;
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
set(gca,'FontSize',18)

saveas(gcf, 'Convergence_Torus', 'epsc');

! rm *.txt *.end *.h5 Run_* .*








