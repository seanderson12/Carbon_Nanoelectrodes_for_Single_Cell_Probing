clf;
nuc=1;                                                                                 %nuc=1 for Nuclear Membrane Inclusion 0 for Cyt. Only                                                                
    n1=1;
    for dpthend=1e-7:.1e-7:3.5e-6;                                                     %Parameter Sweep 1, Final Pipette Depth in Cell (m)                                                             
            %Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            r_tip=125e-9;                                                              %CNP Tip Radius (m)
            LCNP=4e-6;                                                                 %CNP Length (m)
            Rc=(22e-6);                                                                %Cell Radius (m)
            Rn=(5e-6);                                                                 %Nuclear Radius (m)
            RS=10e6;                                                                   %Series Resistance (Ohm), Exp. Approximation
            RSc=1e5;                                                                   %Intracellular Series Resistance (Ohm), Exp. Approximation
            Rct=2e9;                                                                   %Charge Transfer Resistance Outside Cell (Ohm), Exp Approximation                                                                 
            RM=.5e9;                                                                   %Membrane Resistance (Ohm), Lit. Approximation
            AR=(Rn^2)/(Rc)^2;                                                          %Membrane Area ratio
            RM2=RM/AR;                                                                 %Nuclear Membrane Resistance (Ohm)
            co=0.16;                                                                   %Molar concentration outside of the cell, Exp. Data
            ci=0.155;                                                                  %Molar concentration inside the cell, Lit. Data            
            
            %CNP Area Calculation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            At=(pi*(LCNP+r_tip/tan(1.6*pi()/180))*tan(1.6*pi/180)...
                *(LCNP+r_tip/tan(1.6*pi()/180))/cos(1.6*pi/180))...
                -pi*r_tip^2/sin(1.6*pi/180);                                           %Total CNP Area 1.6 Degree Cone Angle (Exp Meas.)
            dpth=linspace(1e-7,dpthend,10);                                            %Depth inside Cell                         
            Ai=pi.*(dpth+r_tip/tan(1.6*pi()/180)).*tan(1.6*pi/180)...
                .*(dpth+r_tip/tan(1.6*pi()/180))./cos(1.6*pi/180)...
                -pi*r_tip^2/sin(1.6*pi/180);                                           %Area Inside Cell           
            Ao=At-Ai;                                                                  %Area Outside Cell                                                      
            s=size(dpth);                                                              %Size of dpth Vector  
            
            %Capacitance Components%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Am=4*pi()*(Rc)^2;                                                          %Memrane Area Assuming Spherical Cell of Radius Rc
            CM=(1e-6)*Am*(10^(4));                                                     %Membrane Capacitance (1uF/cm^2)
            Csa=(64.74);                                                               %Stern Capacitance (uF/cm^2), Theory 450pm Hydrated Radius Sodium, and Rel Permittivity of 78.5 for Media
            Cis=(Csa*10^(-6))*Ai*(10^4);                                               %Stern Layer Capacitance Inner 
            Cos=(Csa*10^(-6))*Ao*(10^4);                                               %Stern Layer Capacitance Outer 
            Ci=((228*(ci^(1/2))*(10^-2).*Ai*cosh(19.5*0.01))...
                .^(-1)+Cis.^(-1)).^(-1);                                               %EDL Capacitance Inside Cell (F)
            Co=((228*(co^(1/2))*(10^-2).*Ao*cosh(19.5*0.01))...
                .^(-1)+Cos.^(-1)).^(-1);                                               %EDL Capacitance Outside Cell (F)
            CM2=(1e-6)*AR*Am*(10^4);
            C_Stray=0;                                                                 %Stray Capacitance from the System (F)                                                     
            C_Inner=0;                                                                 %Inner Capacitance of Pipette (F)
            Rcti=Rct*At./Ai;
            Rcto=Rct*At./Ao;
            
            %TOTAL IMPEDANCE CALC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            w=2*pi*1000;                                                               %Angular Frequency at 1000Hz
            Zo=((1./Rcto)+(1i*w.*Co)).^(-1);                                           %Extracellular Electrode Impedance
            Zi=((1./Rcti)+(1i*w.*(Ci+C_Inner))).^(-1)+RSc;                             %Intracellular Electrode Impedance
            ZM=((1/RM)+(1i*w.*CM)).^(-1);                                              %Membrane Impedance
            if nuc==1                                                                  %Possibility for Inclusion of nuclear membrane
                ZM2=((1/RM2)+(1i*w.*CM2)).^(-1);
                RSc=2*RSc;
            else
                ZM2=0;
            end
            ZS=RS;                                                                     %Series Impedance
            Znet=ZS+(((Zi+ZM+ZM2).^(-1)+(Zo).^(-1)).^(-1));                            %Net Impedance
            Re=real(Znet);                                                             %Real Part of Impedance
            Im=imag(Znet);                                                             %Imaginary Part of Impedance
            RDC=RS+((1./Rcto)+(1./(Rcti+RM+RSc))).^(-1);                               %DC Resistance
            A=real(1./Znet);                                                           %Real Part of Admittance
            B=imag(1./Znet);                                                           %Imaginary Part of Admittance
            b=1./RDC;                                                                  %DC Conductance
            CM_EPC10=(1./(w.*B)).*((A.^2+B.^2-A.*b).^2)./((A-b).^2+B.^2);              %Sine+DC Algorithm Capacitance (algorithm used by amplifier software in experiment)
            DRE(n1)=Re(end)-Re(1);                                                     %Change in Real Part of Impedance
            DIM(n1)=Im(end)-Im(1);                                                     %Change in Imag Part of Impedance
            DC(n1)=CM_EPC10(end)-CM_EPC10(1);                                          %Change in Measured Capacitance
        n1=n1+1;                                                                       %Incrementing index
    end

%%%%%%%%%%PLOTTING DATA%%%%%%%%%%%%%%%%
dpthend=1e-7:.1e-7:3.5e-6;
DC=DC*10^15;                                                                           %Converting Cap. from F-->fF
DRE=DRE/1000;                                                                          %Converting Imp. from Ohm-->kOhm
DIM=DIM/1000;                                                                          %Converting Imp. from Ohm-->kOhm
h5=plot((dpthend*10^6),DC,'-.k','LineWidth',2);
hold on
h6=plot((dpthend*10^6),DRE,'-k','LineWidth',2);
set(gca,'FontSize',15)
legend('Capacitance Change (fF)','Re(Z) Change (kOhm)',...
    'Location','NorthWest');
xlabel('Depth in Cell (um)');
ylabel('Ceq (fF) or Re(Z) (kOhm) Change');
