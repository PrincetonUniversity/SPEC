function plot_spec_Jvolume(filename, component, ns, theta, zeta, newfig)


% Constant definition
mu0 = 4*pi*1E-7;
epsilon = 1E-5;

% Data loading
fdata = read_spec_field(filename);      % Read data
Nvol = fdata.Nvol;                      % Total number of volumes
sarr = linspace(-1, 1, ns);

% Allocate memory
Bcontrav = cell(1, Nvol);
Bcov = cell(1, Nvol);

% Get magnetic field
for ivol=1:Nvol
    if ivol==1
        sarr(1)=-1+epsilon;
    else
        sarr(1)=-1;
    end
    temp = get_spec_magfield(fdata, ivol, sarr, theta, zeta);
    Bcontrav{ivol} = zeros(3, length(temp{1}));
    
    for ii=1:3 % transform the output of get_spec_magfield_cyl in an array
        Bcontrav{ivol}(ii,:) = temp{ii};
    end
    
    Bcov{ivol} = contra2cov(filename, ivol, sarr, Bcontrav{ivol}, ns, ...
                            theta, zeta, 0);
end

% Data processing

% First, get the current in each volume
nr = Nvol*ns;                       % Number of points in r coordinate
j_parallel_contrav = zeros(3, nr);  % Allocate memory
j_parallel_covar = zeros(3, nr);    % Allocate memory
r = zeros(1, nr);                   % Allocate memory
r_kam = zeros(1, Nvol);             % Allocate memory

iimin = 1;
iimax = ns;
for ivol=1:Nvol

    % Construct radial coordinate
    if ivol==1
        sarr(1)=-1+epsilon;
        sbar = (sarr+1) / 2;
        rmax = get_spec_radius(filename, theta, zeta, ivol);
        r(iimin:iimax) = rmax .* sqrt(sbar);
    else
        sarr(1)=-1;
        sbar = (sarr+1) / 2;
        rmin = get_spec_radius(filename, theta, zeta, ivol-1);
        rmax = get_spec_radius(filename, theta, zeta, ivol);
        r(iimin:iimax) = rmin + (rmax-rmin) * sbar;
    end
    
    r_kam(ivol) = rmax;

    % And compute the parallel current
    j_parallel_contrav(:,iimin:iimax) = fdata.mu(ivol)/ mu0 * Bcontrav{ivol};
    j_parallel_covar(:,iimin:iimax) = contra2cov(filename, ivol, sarr,...
                          j_parallel_contrav(:,iimin:iimax), ns, theta, zeta, 1);

    % Change indices for next volume
    iimin = iimin + ns;
    iimax = iimax + ns;
end




% some plots

switch newfig
    case 0
        hold on
    case 1
        figure
        hold on
    case 2
        hold off
end

switch component
    case 'psi'
        plot(r, j_parallel_covar(1,:), 'LineWidth', 1.5, 'DisplayName', '$j_{\parallel,\psi}$')
        hold on;
    case 'theta'
        plot(r, j_parallel_covar(2,:), 'LineWidth', 1.5, 'DisplayName', '$j_{\parallel,\theta}$')
        hold on;
    case 'phi'
        plot(r, j_parallel_covar(3,:), 'LineWidth', 1.5, 'DisplayName', '$j_{\parallel,\phi}$')
        hold on;
    case 'all'
        plot(r, j_parallel_covar(1,:), 'LineWidth', 1.5, 'DisplayName', '$j_{\parallel,\psi}$')
        hold on;
        plot(r, j_parallel_covar(2,:), 'LineWidth', 1.5, 'DisplayName', '$j_{\parallel,\theta}$')
        plot(r, j_parallel_covar(3,:), 'LineWidth', 1.5, 'DisplayName', '$j_{\parallel,\phi}$')
end
        
leg = legend('Location','northwest');
ylab = ylabel('$\mathbf{J}_\mathcal{V}$');
xlab = xlabel('r');
set(gca, 'FontSize', 14)
set(leg,'Interpreter','latex');
set(xlab,'Interpreter','latex');
set(ylab,'Interpreter','latex');
grid on;
