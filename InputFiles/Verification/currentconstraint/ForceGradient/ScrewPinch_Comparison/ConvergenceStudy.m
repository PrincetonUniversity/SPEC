
% Case Nvol=2

fname = 'ScrewPinch_Nvol3.sp';
Lrad  = 6:16;

system(['../../../../../dspec ', fname])
outname = [fname, '.h5'];
data    = read_spec(outname);
Nvol   = data.input.physics.Nvol;


ErrHessian = zeros(1,length(Lrad));

hessians = cell(1, length(Lrad));
for ii=1:length(Lrad)
    disp(newline)
    disp(['Run ', num2str(ii), ' / ', num2str(length(Lrad))])
    disp(newline)
    
    InputFile = ['Run_', num2str(ii), '.sp'];
    
    change_spec_inputfile(fname, InputFile, 'Lrad', ones(1,Nvol) * Lrad(ii), ...
                          'Lfindzero', 0, 'Lcheck', 0, 'LHmatrix', 'T', 'Mregular', 4);  
    system(['../../../../../dspec ', InputFile])
    
    hessians{ii} = read_spec_hessian([InputFile, '.h5']);
    
    data = read_spec([InputFile, '.h5']);
    Nvol  = data.input.physics.Nvol;
    tflux = data.output.tflux;
    pflux = data.output.pflux;
    mu    = data.output.mu;
    R     = data.output.Rbc;
    ns    = 2560000;             

    tflux(2:end) = diff(tflux);
    pflux(2:end) = diff(pflux);

    hessian_analytic = get_hessian_screwpinch(Nvol, tflux, pflux, mu, R, ns);
    
    ErrHessian(ii) = max(max(abs(hessians{ii} - hessian_analytic)));
    
end


figure
plot(Lrad, log10(ErrHessian), '*', 'MarkerSize', 8, 'LineWidth', 2)
xlabel('log_{10}(Lrad)')
ylabel('log_{10}(\Delta H)')
legend(['Screw pinch, N_{vol}=',num2str(Nvol)])
set(gca, 'FontSize', 18)
set(gcf, 'color', 'w')
grid on


saveas(gcf,'Convergence_ScrewPinch','epsc')


! rm Run_*
! rm ScrewPinch_Nvol3.sp.*
