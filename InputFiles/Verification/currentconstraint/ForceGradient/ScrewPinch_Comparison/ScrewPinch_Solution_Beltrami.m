function struc = ScrewPinch_Solution_Beltrami(Rmin, Rmax, ns, pflux, tflux, mu, vvol, innout, ivol)

%innout: B_in or Bout
%vvol  : field in volume vvol
%ivol  : label perturbed interface


rarr = linspace(Rmin, Rmax, ns);

if innout==0
    R = Rmin;
else
    R = Rmax;
end

besselj0 = besselj(0,mu*R);
besselj1 = besselj(1,mu*R);
besselj2 = besselj(2,mu*R);

dbesselj0 = -besselj1;
dbesselj1 = 0.5 * (besselj0 - besselj2);



if vvol==1
    if (innout==0)
        error('Incompatible Lsingularity with innout iocons')
    end
    
    
    Jbar = trapz(rarr, rarr.*besselj(0,mu*rarr));
    %syms r
    %J(r) = r * besselj(0, mu* r);
    %Jbar = vpaintegral(J, Rmin, Rmax);
    dJbar = R*besselj0;
    c1 = tflux(1) / (2*pi*Jbar);
    dc1 = -c1 * dJbar / Jbar;
    dc2 = 0;
    c2 = 0;
    
    B  = c1 * sqrt(besselj0^2 + besselj1^2);                           % |B| at interface
    Bt = c1 * R * besselj(1, mu*R);                                        % Poloidal component of B at interface
    Bz = c1 * R * besselj(0, mu*R);                                        % Toroidal component of B at interface
    BdBdc1 = B*B / c1;                                                     % Derivative of |B| w.r.t c1
    BdBdc2 = 0;                                                            % Derivative of |B| w.r.t c2
    BdBdpf = 0;                                                            % |B| times the derivative of |B| w.r.t pflux
    dBtdpf= 0;                                                             % Derivative of Bt w.r.t pflux
    dBzdpf= 0;                                                             % Derivative of Bz w.r.t pflux
    sg    = R^2/4;
    
    switch ivol
        case 1


            BdBdR  = mu * c1^2 * (besselj0 * dbesselj0 + besselj1 * dbesselj1)...% |B| times the derivative of |B| w.r.t R
                   + c1*dc1 * (besselj1^2 + besselj0^2);

            dBtdR =  c1     *  besselj1        ...                                 % Derivative of Bt w.r.t R
                  +  c1 * R * dbesselj1 * mu...
                  + dc1 * R *  besselj1;

            dBzdR =  c1     *  besselj0        ...                                 % Derivative of Bz w.r.t R
                  +  c1 * R * dbesselj0 * mu...
                  + dc1 * R *  besselj0;
            Rs = R/4;
            
            sgBtup = R/4 * c1 * besselj1;
            dsgBtupdR = c1*besselj1/4 + R/4 * dc1 * besselj1 + R/4 * c1 * mu * dbesselj1;
            
            gtt_over_sg = R/Rs;
            dgtt_over_sgdR=0;
            
            dsgdR = R/2;
            
            
            dBtupdR = -1/R^2 * c1 * besselj1 + 1/R * dc1 * besselj1 + 1/R * c1 * mu * dbesselj1;
            
        otherwise
            
            BdBdR = 0;
            dBtdR = 0;
            dBzdR = 0;
            sgBtup = 0;
            dsgBtupdR = 0;
            gtt_over_sg=0;
            dgtt_over_sgdR=0;
            dsgdR  = 0;
            
            dBtupdR = 0;
    end
    
    
    
    
    
else
    
    bessely0 = bessely(0,mu*R);
    bessely1 = bessely(1,mu*R);
    bessely2 = bessely(2,mu*R);

    dbessely0 = -bessely1;
    dbessely1 = 0.5 * (bessely0 - bessely2);
    
    Jbar = trapz(rarr, rarr.*besselj(0,mu*rarr));
    Ybar = trapz(rarr, rarr.*bessely(0,mu*rarr));

    Jtil = trapz(rarr,       besselj(1,mu*rarr));
    Ytil = trapz(rarr,       bessely(1,mu*rarr));
    
    %syms r;
    %J0(r) = r*besselj(0,mu*r);
    %Y0(r) = r*bessely(0,mu*r);
    %J1(r) =   besselj(1,mu*r);
    %Y1(r) =   bessely(1,mu*r);
    
    %Jbar = vpaintegral(J0, Rmin, Rmax);
    %Ybar = vpaintegral(Y0, Rmin, Rmax);
    %Jtil = vpaintegral(J1, Rmin, Rmax);
    %Ytil = vpaintegral(Y1, Rmin, Rmax);
    
    
    
    D = Jbar * Ytil - Jtil * Ybar;
    c1  = (Ytil*tflux - Ybar*pflux) / (2.0*pi*D);
    c2  = (Jbar*pflux - Jtil*tflux) / (2.0*pi*D);
    
    sg = R * (Rmax-Rmin)/2.0;
    
    switch ivol
        case vvol-1
            dJbar = -besselj(0,mu*Rmin) * Rmin;
            dYbar = -bessely(0,mu*Rmin) * Rmin;
            dJtil = -besselj(1,mu*Rmin);
            dYtil = -bessely(1,mu*Rmin);
        case vvol
            dJbar = besselj(0,mu*Rmax) * Rmax;
            dYbar = bessely(0,mu*Rmax) * Rmax;
            dJtil = besselj(1,mu*Rmax);
            dYtil = bessely(1,mu*Rmax);            
        otherwise
            dJbar = 0;
            dYbar = 0;
            dJtil = 0;
            dYtil = 0;           
            
    end
    c  = [c1;c2];
    dJ = [[dJbar, dYbar]; [dJtil, dYtil]];
    J  = [[Jbar, Ybar]; [Jtil, Ytil]];

    dc = -J\(dJ*c);
    dc1 = dc(1);
    dc2 = dc(2);
    
    
    Bt = R * (c1 * besselj1 + c2 * bessely1);                              % Poloidal component of B at interface
    Bz = R * (c2 * besselj0 + c2 * bessely0);                              % Toroidal component of B at interface

    if innout==0 % evaluate at inner interface (Rmin)
        switch ivol
            case vvol-1 % Take derivative w.r.t Rmin
                dBtdR   = c1 * besselj1 + c2 * bessely1 + mu * R * (c1 * dbesselj1 + c2 * dbessely1)...% Derivative of Bt w.r.t R
                        + R * dc1 * besselj1 + R * dc2 * bessely1;
                dBzdR   = c1 * besselj0 + c2 * bessely0 + mu * R * (c1 * dbesselj0 + c2 * dbessely0)...% Derivative of Bz w.r.t R
                        + R * dc1 * besselj0 + R * dc2 * bessely0;
                BdBdR   = c1^2 * mu * ( besselj0 * dbesselj0 + besselj1 * dbesselj1)...% |B| times the derivative of |B| w.r.t R
                        + c2^2 * mu * ( bessely0 * dbessely0 + bessely1 * dbessely1)...
                        + c1*c2* mu * (dbesselj0 * bessely0  + besselj0 * dbessely0 ...
                                     + dbesselj1 * bessely1  + besselj1 * dbessely1)...
                        + c1*dc1*(besselj0^2 + besselj1^2) + c2*dc2 * (bessely0^2 + bessely1^2)...
                        + (dc1/c1 + dc2/c2)*c1*c2*(besselj0*bessely0 + besselj1*bessely1);
                dsgBtupdR = -0.5 * (c1 * besselj1 + c2 * bessely1)...
                          + (Rmax-Rmin)/2.0 * (dc1 * besselj1 + dc2 * bessely1)...
                          + (Rmax-Rmin)/2.0 * mu * (c1 * dbesselj1 + c2 * dbessely1);
                dgtt_over_sgdR = 2.0/(Rmax-Rmin) + 2*Rmin / (Rmax-Rmin)^2;
                
                dsgdR   = (Rmax-Rmin)/2.0 - Rmin/2.0;
                
                dBtupdR = -(c1*besselj1+c2*bessely1)/R^2 + 1/R*(dc1*besselj1 + dc2*bessely1)...
                        + mu/R*(c1*dbesselj1 + c2*dbessely1);
                
            case vvol % Take derivative w.r.t Rmax
                dBtdR   = Rmin * dc1 * besselj1 + Rmin * dc2 * bessely1;
                dBzdR   = Rmin * dc1 * besselj0 + Rmin * dc2 * bessely0;
                BdBdR   = c1*dc1*(besselj0^2 + besselj1^2) ...
                        + c2*dc2*(bessely0^2 + bessely1^2)   ...
                        + (dc1/c1 + dc2/c2)*c1*c2*(besselj0*bessely0 + besselj1*bessely1);
                dsgBtupdR = 0.5 * (c1 * besselj1 + c2 * bessely1) + (Rmax-Rmin)/2.0 * (dc1*besselj1 + dc2 * bessely1);
                dgtt_over_sgdR = - 2*Rmin / (Rmax-Rmin)^2;
                
                dsgdR   = Rmin / 2.0;
                
                dBtupdR = 1/Rmin*(dc1*besselj1 + dc2*bessely1);
            otherwise
                dBtdR   = 0;
                dBzdR   = 0;
                BdBdR   = 0;
                dsgBtupdR = 0;
                dgtt_over_sgdR = 0;
                dsgdR   = 0;
                dBtupdR = 0;
        end
    else % evaluate at outer interface (Rmax)
       switch ivol
           case vvol-1 % Take derivative w.r.t Rmin
                dBtdR   = Rmax * dc1 * besselj1 + Rmax * dc2 * bessely1;
                dBzdR   = Rmax * dc1 * besselj0 + Rmax * dc2 * bessely0;
                BdBdR   = c1*dc1*(besselj0^2 + besselj1^2) + c2*dc2*(bessely0^2 + bessely1^2)...
                        + (dc1/c1 + dc2/c2)*c1*c2*(besselj0*bessely0 + besselj1*bessely1);
                dsgBtupdR = -0.5 * (c1 * besselj1 + c2 * bessely1) + (Rmax-Rmin)/2.0 * (dc1*besselj1 + dc2 * bessely1);
                dgtt_over_sgdR = 2*Rmax / (Rmax-Rmin)^2;
                
                dsgdR   = -Rmax / 2.0;
                dBtupdR = 1/Rmax*(dc1*besselj1 + dc2*bessely1);
                
           case vvol % Take derivative w.r.t Rmax
                dBtdR   = c1 * besselj1 + c2 * bessely1 + mu * R * (c1 * dbesselj1 + c2 * dbessely1)...% Derivative of Bt w.r.t R
                        + R * dc1 * besselj1 + R * dc2 * bessely1;
                dBzdR   = c1 * besselj0 + c2 * bessely0 + mu * R * (c1 * dbesselj0 + c2 * dbessely0)...% Derivative of Bz w.r.t R
                        + R * dc1 * besselj0 + R * dc2 * bessely0;
                BdBdR   = c1^2 * mu * ( besselj0 * dbesselj0 + besselj1 * dbesselj1)...% |B| times the derivative of |B| w.r.t R (at c_i=cst)
                        + c2^2 * mu * ( bessely0 * dbessely0 + bessely1 * dbessely1)...
                        + c1*c2* mu * (dbesselj0 * bessely0  + besselj0 * dbessely0 ...
                                     + dbesselj1 * bessely1  + besselj1 * dbessely1)...
                        + c1*dc1*(besselj0^2 + besselj1^2) + c2*dc2*(bessely0^2 + bessely1^2)...
                        + (dc1/c1 + dc2/c2)*c1*c2*(besselj0*bessely0 + besselj1*bessely1);
                dsgBtupdR = 0.5 * (c1 * besselj1 + c2 * bessely1)...
                          + (Rmax-Rmin)/2.0 * (dc1 * besselj1 + dc2 * bessely1)...
                          + (Rmax-Rmin)/2.0 * mu * (c1 * dbesselj1 + c2 * dbessely1);
                dgtt_over_sgdR = 2/(Rmax-Rmin) - 2*Rmax / (Rmax-Rmin)^2;
                
                dsgdR   = (Rmax-Rmin)/2.0 + Rmax/2.0;
                
                dBtupdR = -(c1*besselj1+c2*bessely1)/R^2 + 1/R*(dc1*besselj1 + dc2*bessely1)...
                        + mu/R*(c1*dbesselj1 + c2*dbessely1);
           otherwise
                dBtdR   = 0;
                dBzdR   = 0;
                BdBdR   = 0;
                dsgBtupdR = 0;
                dgtt_over_sgdR = 0;
                dsgdR   = 0;
                
                dBtupdR = 0;
       end        
    end
    
    dBtdpf = R * (-Ybar * besselj1 + Jbar * bessely1) / (2*pi*D);          % Derivative of Bt w.r.t pflux
    dBzdpf = R * (-Ybar * besselj0 + Jbar * bessely0) / (2*pi*D);          % Derivative of Bz w.r.t pflux
    
    B       = sqrt(c1^2 * ( besselj0^2 + besselj1^2) ...                   % |B| at interface
            + c2^2 * (bessely0^2 + bessely1^2)       ...
            + 2*c1*c2*(besselj0*bessely0 + besselj1*bessely1));
    
    BdBdc1 = c1 * (besselj0^2 + besselj1^2)...                             % |B| times the derivative of |B| w.r.t c1
           + c2 * (besselj0*bessely0 + besselj1*bessely1);
    BdBdc2 = c2 * (bessely0^2 + bessely1^2)...                             % |B| times the derivative of |B| w.r.t c2
           + c1 * (besselj0*bessely0 + besselj1*bessely1);

    dc1dpf = - Ybar / (2*pi*D);
    dc2dpf =   Jbar / (2*pi*D);

    BdBdpf = BdBdc1 * dc1dpf + BdBdc2 * dc2dpf;                            % |B| times the derivative of |B| w.r.t pflux
    
    
    
    
    
    
    gtt_over_sg = 2*R / (Rmax-Rmin);    
    sgBtup = (Rmax-Rmin)/2.0 * Bt / R;
    
    
    
    %test
    tmp1 = dgtt_over_sgdR * sgBtup + gtt_over_sg * dsgBtupdR;
    tmp2 = dBtdR;
    
    %disp(['Direct: ', num2str(tmp2,16), ' or via contravariant: ', num2str(tmp1,16)])
    
    %disp(newline)
    
    
end



struc.B         = B;
struc.BdBdR     = BdBdR;
struc.BdBdc1    = BdBdc1;
struc.BdBdc2    = BdBdc2;
struc.BdBdpf    = BdBdpf;
struc.Bt        = Bt;
struc.Bz        = Bz;
struc.dBtdR     = dBtdR;
struc.dBtdpf    = dBtdpf;
struc.dBzdR     = dBzdR;
struc.dBzdpf    = dBzdpf;
struc.c1        = c1;
struc.c2        = c2;
struc.dc1       = dc1;
struc.dc2       = dc2;
struc.gtt_over_sg = gtt_over_sg;
struc.dgtt_over_sgdR = dgtt_over_sgdR;
struc.sgBtup = sgBtup;
struc.dsgBtupdR = dsgBtupdR;
struc.sg        = sg;
struc.dsgdR     = dsgdR;
struc.dBtupdR   = dBtupdR;

end 
