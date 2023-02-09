%% SCRIPT SCANNING R10 and R11 FROM SPEC AND EXTRACTING MAGNETIC SHEAR %%
%-----------------------------------------%
% Created by Salomon Guinchard (04/26/22) %
% Last modified (10/10/22)                %
%-----------------------------------------%

%% Paths to get_spec_mat_iota_torsion and other SPEC analysis routines (change path according to output folder)
addpath(genpath('/path_to_OutputSPEC'))
addpath(genpath('/path_to_/OutputBOZ'))
fs=18; lw=1.5; % Fontsize, linewidth


%% PARAMETERS %%
nscan = 20;
f = strings(nscan,nscan);
nn = length(f(:,1));
parfor i=1:nn
    for j=1:nn
        f(i,j) = ['File', num2str(i), '_', num2str(j) '.sp.h5'];
    end 
end
filenames = reshape(f,1,nn^2);
ll = length(filenames);
d0 = read_spec(char(filenames(1)));

l = d0.transform.fiota(:,2);
avg_shear   = zeros(1,ll) ;
shift = 5;

disp 'Done initialising...'
%% DISCRETISATION OF SPACE %%
scan__11 = linspace(-0.7,0.7,nscan); % values between which the mode 11 will be scanned
scan__10 = linspace(-3,3,nscan);     % same for 10 mode (01 - see SPEC input format)

%% RUN %%
tic % set timer on
parfor ii=1:ll
    try 
    d = read_spec(char(filenames(ii)));
        phi = 0; 

        Nfp  = double(d.input.physics.Nfp);   % Number of field periods
        Ntor = double(d.input.physics.Ntor);  % Number of toroidal planes  
        Nppts = d.input.diagnostics.nPpts;    % # points Poincare plots

        R0n = double(d.output.Rbc(1:Ntor+1,1));
        Z0n = double(d.output.Zbs(1:Ntor+1,1));

        % Define coordinate axis position
        Raxis = double(0);
        Zaxis = double(0);
        for i=0:Ntor
            Raxis = Raxis + R0n(i+1) * cos(i * phi);
            Zaxis = Zaxis + Z0n(i+1) * sin(i * phi);
        end

        Pos_Axis = [Raxis Zaxis]; 

        tt = length(d.poincare.R(:,1,1));
            jj=1;
            X =  reshape(d.poincare.R(jj,1,:),[1 ,Nppts]);
            Y =  reshape(d.poincare.Z(jj,1,:),[1 ,Nppts]);
            Xbar = mean(X);
            Ybar = mean(Y);
            theta = atan2(Y-Ybar , X-Xbar);
            theta = mod(theta, 2*pi);

            [theta, ind] = sort(theta);

            X = X(ind);
            Y = Y(ind);
            p  = polyshape(X,Y);
            %clc
            is = isinterior(p,Pos_Axis);

        while is == 0
             jj = jj+1;

        X =  reshape(d.poincare.R(jj,1,:),[1 ,Nppts]);
        Y =  reshape(d.poincare.Z(jj,1,:),[1 ,Nppts]);

        Xbar = mean(X);
        Ybar = mean(Y);
        theta = atan2(Y-Ybar , X-Xbar);
        theta = mod(theta, 2*pi);

        [theta, ind] = sort(theta);

        X = X(ind);
        Y = Y(ind);
        p  = polyshape(X,Y);
        %clc
        is = isinterior(p,Pos_Axis);
        
        end
        n_surf = jj;
        n = length(l) - (n_surf+shift-1);
        
        mat_iota    = zeros(1,n);
        mat_r_coord = zeros(1,n);
        mat_s_coord = zeros(1,n);
        derivatives = zeros(1,n);
        shear       = zeros(1,n);
        coeff       = zeros(1,n);

    out = extract_shear(d,n,shift, n_surf);
    mat_iota          = out.mat_iota;
    mat_r_coord       = out.mat_r_coord;
    mat_s_coord       = out.mat_s_coord;
    scan_11(ii)       = out.scan_11;
    scan_10(ii)       = out.scan_10;
    derivatives       = out.derivatives;
    avg_shear(ii)     = out.avg_shear;
    catch 
    warning('Problem extracting shear.  Assigning a value of -1 to d and all parameters to NaN');
    d = -1;
    n =  1;
    mat_iota          = NaN(1,n);
    mat_r_coord       = NaN(1,n);
    mat_s_coord       = NaN(1,n);
    scan_11(ii)       = NaN;
    scan_10(ii)       = NaN;
    derivatives       = NaN(1,n);
    avg_shear(ii)     = NaN;
    end

end

disp 'Done running...'
toc % set timer off
%% Plot last polyshape and coordinate axis %%
    figure('color','w')
    pg=plot(p);
    pg.FaceColor = ([1 1 1]);
    pg.EdgeColor = ([1 0 0]);
    pg.LineWidth =   3;
    hold on 
    plot(Pos_Axis(1), Pos_Axis(2), 'k+', 'linewidth' , 3, 'MarkerSize',8)
    set (gca, 'fontsize', 16)
    xlabel('R', 'interpreter', 'latex', 'FontSize', 24)
    ylabel('Z', 'interpreter', 'latex','FontSize', 24)

%% Shear %%

shear_mat = reshape(avg_shear,[nn,nn])';
r10_tilde = reshape(scan_10,  [nn,nn]);
r11_tilde = reshape(scan_11,  [nn,nn]);
parfor i=1:nn
     r11(i) = r11_tilde(1,i);
     r10(i) = r10_tilde(i,1);
     
    for j=1:nn
        S(i,j) = shear_mat(i,j);
    end 
end 
[R10,R11] = meshgrid(r10,r11);


%% PCOLOR SHEAR %%
figure
pcolor(R11,R10,S)
shading interp
hold on 
contour(R11,R10,S,[0,0], 'r', 'LineWidth', 3)
set (gca, 'fontsize', fs)
colorbar
xlabel('$R_{11}$', 'interpreter', 'latex')
ylabel('$R_{10}$', 'interpreter', 'latex')


%% PLOT PROFIL IOTA ACROSS ZERO SHEAR LEVEL CURVE %%

filename_1 = ['File1', num2str(index_R10_1), '_' num2str(index_R11) '.sp.h5'];
filename_2 = ['File2', num2str(index_R10_2), '_' num2str(index_R11) '.sp.h5'];

shift = 8; % shift in extract shear function. Has to be determined from iota profile

d_bel  = read_spec(filename_1);
d_abov = read_spec(filename_2);

coord_1 = d_bel.transform.fiota(shift:end,1);
coord_2 = d_abov.transform.fiota(shift:end,1);

iota_1 = d_bel.transform.fiota(shift:end,2);
iota_2 = d_abov.transform.fiota(shift:end,2);

figure
plot(coord_1,iota_1,'b+', 'linewidth', 2, 'MarkerSize', 7)
grid on
hold on
plot(coord_2,iota_2,'r+', 'linewidth', 2, 'MarkerSize', 7)
set(gca, 'fontsize', 17)
xlabel('$r$', 'Interpreter','latex', 'FontSize',30)
ylabel('$\iota$', 'Interpreter','latex', 'FontSize',30)
legend('$\iota_{a}$', '$\iota_{b}$', 'Interpreter', 'latex',   'FontSize', 27)


%% PCOLOR & SURF MULTIPLOT %%

figure
subplot(1,2,1)
%----------------------------------------------------------
pcolor(R11,R10, S)
title('Shear($R_{11},R_{10}$)', 'interpreter', 'latex')
shading interp
hold on 
contour(R11,R10,S,[0,0], 'r', 'LineWidth', 3)
hold on
plot(scan11_trunc, scan10_trunc, 'w--', 'linewidth', 1);
hold on
plot(scan11_trunc2, scan10_trunc2, 'w--', 'linewidth', 1);
set (gca, 'fontsize', fs)
colorbar
xlabel('$R_{11}$', 'interpreter', 'latex')
ylabel('$R_{10}$', 'interpreter', 'latex')

subplot(1,2,2)
%---------------------------------------------------------
surf(R11,R10,S)
shading interp
hold on 
contour(R11,R10,S,[0,0], 'r', 'LineWidth', 3)
set (gca, 'fontsize', fs)
colorbar
xlabel('$R_{11}$', 'interpreter', 'latex')
ylabel('$R_{10}$', 'interpreter', 'latex')

clc


%% Extract %%
function out = extract_shear(d, n, shift, n_surf) % can be found in extract_shear.m file too
    
  id = n_surf+shift;
  radial_coord     = (d.transform.fiota(:,1));       % extract radial coordinate
  out.mat_iota     =  d.transform.fiota(id:end,2);   % extract iota
  out.mat_r_coord  = radial_coord(id:end);           % truncates the radial coordinates (remove 5 first terms)
  out.mat_s_coord  = ((out.mat_r_coord + 1)./2);     % change of variable r <--> s
  out.scan_11      = d.output.Rbc(11,2);             % Value of R11
  out.scan_10      = d.output.Rbc(2,2);              % Value of R10
 
  % CENTERED FINITE DIFFERENCES 
    out.derivatives(1) = (out.mat_iota(2) - out.mat_iota(1))   / (out.mat_s_coord(2) - out.mat_s_coord(1));
    out.derivatives(n) = (out.mat_iota(n) - out.mat_iota(n-1)) / (out.mat_s_coord(n) - out.mat_s_coord(n-1));
    
 for j=2:n-1
    out.derivatives(j) = (out.mat_iota(j+1) - out.mat_iota(j-1)) / (out.mat_s_coord(j+1) - out.mat_s_coord(j-1));
 end
 
 out.coeff      =    out.mat_s_coord ./ (out.mat_iota) ;
 out.shear      =    (out.coeff)'    .* out.derivatives;
 out.avg_shear  =    mean(out.shear);                     % Avg shear 
end