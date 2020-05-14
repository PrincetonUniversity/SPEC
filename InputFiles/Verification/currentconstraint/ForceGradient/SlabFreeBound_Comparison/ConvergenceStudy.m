

% Case Nvol=2

fname = 'Slab_FreeBound_Nvol1.sp';
Lrad  = 2:16;

system(['../../../../../dspec ', fname])
outname = [fname, '.h5'];
data    = read_spec(outname);

Nvol  = data.input.physics.Nvol;
Mvol  = Nvol + 1;
tflux = data.output.tflux;
pflux = data.output.pflux;
mu    = data.output.mu;
R     = data.output.Rbc;
Isurf = data.output.IPDt;
Icoil = data.input.physics.curpol;

tflux(2:end) = diff(tflux);
pflux(2:end) = diff(pflux);

hessian_analytic = get_hessian_slab(Nvol, tflux, pflux, mu, R, Isurf, Icoil, false);

hessians = cell(1, length(Lrad));
for ii=1:length(Lrad)
    disp(newline)
    disp(['Run ', num2str(ii), ' / ', num2str(length(Lrad))])
    disp(newline)
    
    InputFile = ['Run_', num2str(ii), '.sp'];
    
    change_spec_inputfile(fname, InputFile, 'Lrad', ones(1,Mvol) * Lrad(ii), ...
                          'Lfindzero', 0, 'Lcheck', 0, 'LHmatrix', 'T', 'Mregular', 4);
                                        
    system(['../../../../../dspec ', InputFile])
    
    hessians{ii} = read_spec_hessian([InputFile, '.h5']);
    
end

ErrHessian = zeros(1,length(Lrad));
for ii = 1:length(Lrad)
    
    ErrHessian(ii) = max(max(abs(hessians{ii} - hessian_analytic)));
    
end


figure('Position', [200 200 900 700])
plot(Lrad, log10(ErrHessian), '*', 'MarkerSize', 12, 'LineWidth', 2)
hold on;
plot(xlim, [-15, -15], 'k--', 'LineWidth', 3)
annotation('textbox',[.2 .15 .4 .1],'String','Machine Precision','FitBoxToText','on','FontSize',16, 'EdgeColor', 'w');
xlabel('Lrad')
ylabel('log_{10}(\nabla F_{SPEC} - \nabla F_{analytic})')
set(gca, 'FontSize', 18)
set(gcf, 'color', 'w')
grid on

saveas(gcf, 'Convergence_Slab', 'epsc');

! rm Run_*
! rm *.h5
! rm *.end
! rm .Run*
! rm .Slab*
