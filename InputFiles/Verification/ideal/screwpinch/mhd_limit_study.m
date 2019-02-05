% Run many SPEC simulations with increasing number of volumes in
% cylindrical coordinates, and compare it to the general screw pinch MHD
% solution
%
% Written by A.Baillod (2019)


newsim = true;
% Some constants (do not touch!)
mu0 = 4*pi*1E-7;
nptr = 0;
L = 2*pi;
a = 1;

% =========================================================================
%                              INITIALIZATION
%                              --------------

Nvol = [2, 8, 32]; % Modify to chose the number of volumes used in each run
lrad = 12;

pmax = 0.05/mu0;   % Maximum pressure
dpmax = -5E5;      % Maximum pressure gradient
Bz0 = 1;

% Definition of p(r), iota(r) (Input of screwpinch solver)
lambda = -dpmax / pmax;        
r_mhd = linspace(0, 1, 1E3);
p_mhd = pmax * (1 + tanh(lambda*(a/2-r_mhd)));
iota = 1.0 ./ (1 + r_mhd);


Nvollegend = cell(1, length(Nvol)+1);
Nvollegend{1} = 'MHD';
for i=1:length(Nvol)
   Nvollegend{i+1} = char(['Nvol = ', num2str(Nvol(i))]); 
end


% =========================================================================
%                              SOLVERS
%                              -------

% ScrewPinch solver
[r_out, jp, jz, Bp, Bz] = ScrewPinch_IotaSolver(L, a, r_mhd, p_mhd, ...
                                                iota * 2*pi, Bz0, 0);

% Interpolate solution on r_mhd
jp = interp1(r_out, jp, r_mhd, 'spline'); 
jz = interp1(r_out, jz, r_mhd, 'spline');  
Bp = interp1(r_out, Bp, r_mhd, 'spline'); 
Bz = interp1(r_out, Bz, r_mhd, 'spline');                                           
                                            
% We have now jp(r), jz(r), Bp(r), Bz(r). We compute the flux psi(r)
psit_mhd = tflux(r_mhd(2:end), Bz(2:end), r_mhd(2:end), 1);
psit_norm = psit_mhd ./ max(psit_mhd);                      % Toroidal flux
psit_norm = [0, psit_norm];
psip_mhd = pflux(r_mhd(2:end), Bp(2:end), r_mhd(2:end), 1); % Poloidal flux

% First guess for mu computed from current profiles
mu = mu0*(jz.*Bz + jp.*Bp) ./ (Bp.^2 + Bz.^2);

% Now run all SPEC simulations
files = cell(1, length(Nvol));


if newsim
    for ii=1:length(Nvol)
        files{ii} = run_from_mhd(Nvol(ii), a, lrad, nptr, r_mhd, ...
                    mu0*p_mhd, mu, iota, psit_mhd, psip_mhd, true, 0);
    end
else
    for ii=1:length(Nvol)
       files{ii} = ['G2L1_Nvol', num2str(Nvol(ii)), '.sp.h5'] ;
    end
end


% =========================================================================
%                           DATA PROCESSING
%                           ---------------

% Get field data from SPEC
B = cell(1, length(Nvol));
for iiv1=1:length(Nvol)
    % get canonical B field for each simulation
    B{iiv1} = get_full_field(char(files{iiv1}), r_mhd, 0, 0, 500, a); 
end

% Allocate memory
iiv1 =  zeros(1, length(Nvol));
fdata = cell(1, length(Nvol));
mu_spec = zeros(1, length(Nvol));
Bp_exact = cell(1, length(Nvol));
Bt_exact = cell(1, length(Nvol));
psip_spec = cell(1, length(Nvol));
psit_spec = cell(1, length(Nvol));
r_kam = cell(1, length(Nvol));
psi_p = cell(1, length(Nvol));
psi_p2 = cell(1, length(Nvol));
psi_t = cell(1, length(Nvol));
psi_t2 = cell(1, length(Nvol));

for isim=1:length(Nvol)
    % Load SPEC output
    fdata{isim} = read_spec_field(char(files{isim}));
    
    % Allocate memory
    r_kam{isim} = zeros(1, Nvol(isim));
    psi_p{isim} = zeros(1, Nvol(isim));
    psi_t{isim} = zeros(1, Nvol(isim));
    
    % Compute fluxes and position of KAM surfaces
    for j=1:Nvol(isim)
       psi_p{isim}(j) = get_spec_polflux_cyl(fdata{isim}, j, 0, -1, 1, 500, 500); 
       psi_t{isim}(j) = get_spec_torflux_cyl(fdata{isim}, j, 0, -1, 1, 500, 500); 
       r_kam{isim}(j) = get_spec_radius(char(files{isim}), 0, 0, j);
    end

    % Get total flux (integral from 0 to KAM surface)
    psi_p2{isim} = zeros(1, Nvol(isim));
    psi_t2{isim} = zeros(1, Nvol(isim));
    for j=1:Nvol(isim)
       psi_p2{isim}(j) = sum(psi_p{isim}(1:j)) ;
       psi_t2{isim}(j) = sum(psi_t{isim}(1:j)) ;
    end

    % Find index of KAM surfaces in r_mhd
    [temp, iiv1(isim)] = min(abs(r_mhd - r_kam{isim}(1)));
    
    % Get mu from SPEC output
    temp = transpose(fdata{isim}.mu);
    mu_spec(isim) = temp(1);

    % Compute analytical fields in first volume...
    Bp_exact{isim} = Bz0 * besselj(1, mu_spec(isim) * r_mhd);
    Bt_exact{isim} = Bz0 * besselj(0, mu_spec(isim) * r_mhd);
end



% =========================================================================
%                               FIGURES
%                               -------


% Pressure plot
figure
plot(psit_norm, p_mhd * mu0, '-k')
hold on;
for i=1:length(Nvol)
    plot_spec_pressure(char(files{i}), 0)
end
lines = findobj(gcf, 'Type', 'Line');
for i=1:length(lines)
   lines(i).LineWidth = 2 ;
end

xlabel('\psi_t / max(\psi_t)')
ylabel('\mu_0 p')
legend(Nvollegend)


% Rotational transform plot
figure
plot(r_mhd, iota, 'LineWidth', 4)
hold on;
for i=1:length(Nvol)
   idata = read_spec_iota(char(files{i}));
   pdata = read_spec_poincare(char(files{i}));
   plot_spec_iota(idata, pdata, fdata{i}, 'i', 'R', 0);
end
xlabel('r')
ylabel('\iota')
legend(Nvollegend)
set(gca, 'FontSize', 12)

% Field plots
figure
plot(r_mhd, Bp, 'LineWidth', 3)         % ScrewPinch solution
hold on;
for ii=1:length(Nvol)
    plot(r_mhd, B{ii}(2,:), '-.', 'LineWidth', 2) %SPEC solution
end
ylabel('B_\theta')
xlabel('r')
h = legend(Nvollegend);
set(gca, 'FontSize', 12)
set(h, 'Location', 'northwest')

figure
plot(r_mhd, Bz, 'LineWidth', 3)         % ScrewPinch solution
hold on;
for ii=1:length(Nvol)
    plot(r_mhd, B{ii}(3,:), '-.', 'LineWidth', 2) %SPEC solution
end
ylabel('B_\zeta')
xlabel('r')
h = legend(Nvollegend);
set(gca, 'FontSize', 12)
set(h, 'Location', 'northeast')



% First volume field plots. In the first volume, the analytical solution of 
% MRxMHD equations is known (bessel equations) and can be compared with 
% SPEC    

for isim=1:length(Nvol)
    figure  
    set(gca, 'DefaultLineLineWidth', 2)
    hold on;
    yyaxis left;
    plot(r_mhd, Bp_exact{isim}, '-')
    yyaxis right;
    plot(r_mhd, Bt_exact{isim}, '-')
    yyaxis left;
    plot(r_mhd(1:20:iiv1(isim)), B{isim}(2,1:20:iiv1(isim)), '*')

    yyaxis right;
    plot(r_mhd(1:20:iiv1(isim)), B{isim}(3,1:20:iiv1(isim)), '*')
    title('First volume, comparison with analytical solution normalized')
    xlabel('r')
    yyaxis left;
    ylabel('B_p')
    yyaxis right;
    ylabel('B_t')
end

    
    
