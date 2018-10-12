clear all
load DATA_REQD.mat
syms C_rc O_d T_rg dFa dFsc
       
fcc_parameters

F_a   = F_a0 + dFa;      F_sc = F_sc0 + dFsc; 
T_oil = T_oil0 ;         kc1 = kc10 ;

T_ri0 = ((c_po*F_oil + lamda*F_oil*c_pd)*T_oil  +  c_ps*F_rc*T_rg)/(c_po*F_oil  +  lamda*F_oil*c_pd  +  c_ps*F_rc);
gama  = delH_f * F_oil / (T_ri0*(F_rc*c_ps + F_oil*c_po  + lamda*F_oil*c_pd));
phi0  = 1-(m*C_rc);
K0    = k0*exp(-E_f/(R*T_ri0));

T_r   =  T_ri0*(1-((gama*y_f0*K0*phi0*(1-exp(-alpa*tc*COR*z_r)))/(alpa + (K0*phi0*(1-exp(-alpa*tc*COR*z_r))))));
Kr    =  k0*exp(-E_f/(R*T_r));
y_f1  =  y_f0*alpa/(alpa + (Kr*phi0*(1-exp(-alpa*tc*COR))));
T_ri1 =  T_ri0*(1-((gama*y_f0*Kr*phi0*(1-exp(-alpa*tc*COR)))/(alpa + (Kr*phi0*(1-exp(-alpa*tc*COR))))));
y_g1  = (1+R_r)*F_gi*((y_f1^I_gi) - y_f1)/(1-I_gi);

sig   =  1.1 + sig2*(T_rg - 873);
delH  =  -h1 - (h2*(T_rg-960)) + 0.6*(T_rg-960)^2 ; 
k     =  k_com * exp(E_cb*((1/960)-(1/T_rg))/R);

C_cat =  kc1 * sqrt(tc*exp(-E_cf/(R*T_ri1))/(C_rc^N));
C_sc  =  C_rc + C_cat; 

f1=(F_sc*(C_sc-C_rc)/W) - k*O_d*C_rc;
f2=(Ra*(O_in - O_d)/Wa) - (((n+2+((n+4)*sig))*k*O_d*C_rc*W)/(4*M_c*Wa*(1+sig)));
f3=(T_ri1*F_sc/W) + (T_a*F_a*c_pa/(W*c_ps)) - (T_rg*(F_sc*c_ps + F_a*c_pa)/(W*c_ps)) - (delH*k*O_d*C_rc/(c_ps*M_c));
F=[f1,f2,f3];
Jac=jacobian(F,[C_rc,O_d,T_rg]); %Finding jacobian of the non linear model

N=length(Time);
C=[0 0 1]; %Measurement model
X0=[0.02 0.01 50]; %Initial state guess
P0=200*eye(3); %Initial error co variance guess
Q=eye(3); %Process noise co variance
R=0.9; %Measurement noise
%initial state guess
Xpre=X0;
Ppre=P0;
Xestimate=[];
for i=1:N-1
    %Propogation step
    tspan=[Time(i) Time(i+1)];
    u=InputU(i,:);
    dFa=u(1);
    dFsc=u(2);
    [t,y] = ode23(@(t,y) fcc_fn_to_solve_odemodel(t,y,u), tspan, [Xpre(1) Xpre(2) Xpre(3)]);
    l=size(y);
    Xup=y(l(1),:);
    Xupd(i,:)=Xup;
    J=eval(subs(Jac,{C_rc,O_d,T_rg},Xup));
    Pup=J*Ppre*transpose(J)+Q;
    %Calculating Kalman gain
    K=Pup*transpose(C)/(C*Pup*transpose(C)+R);
    %Correction step
    Xcor=transpose(Xup)+K*(Y_measured_case1(i)-C*transpose(Xup))
    Pcor=Pup-K*C*Pup;
    Xest(i,:)=Xcor;
    Xpre=Xcor;
    Ppre=Pcor;

end
hold on;
plot(Xest(:,3));
plot(Xupd(:,3),'g');
plot(Y_measured_case1,'r'),legend('State Estimate','Propogated State','Measurement');
hold off;


