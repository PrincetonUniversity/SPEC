function plot_spec_iotakam(data,iorq,xaxis,newfig)

%
% PLOT_SPEC_IOTAKAM( DATA, IORQ, XAXIS, NEWFIG )
% ==============================================
%
% Plots rotational transform on each side of each ideal interface
%
% INPUT
% -----
%   -data  : produced by calling read_spec(fname) 
%   -iorq  : plot iota('i') or safety factor ('q')
%   -xaxis : 'R' plots R-position of interfaces (at phi=0) as the x-axis
%            'f' plots toroidal flux as the x-axis
%            'r' plots sqrt(toroidal flux) as the x-axis
%   -newfig: opens(=1) or not(=0) a new figure
%
% written by J.Loizu (2017) 

    % Check input
    if ~any(iorq==['i','q'])
        error('InputError: invalid iorq')
    end

    if ~any(xaxis==['R','f','r'])
        error('InputError: invalid xaxis')
    end

    switch newfig
        case 0
            hold on
        case 1
            figure('Color','w')
            hold on
        case 2
            hold off
        otherwise
            error('InputError: invalid newfig')
    end

    Nvol   = data.input.physics.Nvol;
    tflux  = data.input.physics.tflux;
    iota   = data.input.physics.iota;
    oita   = data.input.physics.oita;
    Rmn    = data.output.Rbc;
    freeb  = data.input.physics.Lfreebound;

    R0     = zeros(1,Nvol);

    for l=1:Nvol
        R0(l) = sum(Rmn(:,l+1));
    end

    if(iorq=='i')
        F = iota(2:end);
        G = oita(2:end);
        Flabel='\iota';
    elseif(iorq=='q')
        F = 1./iota(2:end);
        G = 1./oita(2:end);
        Flabel='q';
    end


    switch xaxis
        case 'R'
            plot(R0,F,'r+','MarkerSize',6,'LineWidth',2)
            hold on;
            plot(R0,G,'m+','MarkerSize',6,'LineWidth',2)
            xlabel('R')
            ylabel(Flabel)
        case 'f'    
            plot(tflux(1:end-freeb),F,'r+','MarkerSize',6,'LineWidth',2)
            hold on;
            plot(tflux(1:end-freeb),G,'m+','MarkerSize',6,'LineWidth',2)
            xlabel('\Psi / \Psi_{edge}')
            ylabel(Flabel)
        case 'r'
            plot(sqrt(tflux(1:end-freeb)),F,'r+','MarkerSize',6,'LineWidth',2)
            hold on;
            plot(sqrt(tflux(1:end-freeb)),G,'m+','MarkerSize',6,'LineWidth',2)
            xlabel('(\Psi / \Psi_{edge})^{1/2}')
            ylabel(Flabel)  
    end
    
end
