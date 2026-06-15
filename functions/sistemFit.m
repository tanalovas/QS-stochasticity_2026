function dydt = sistemFit(~,y,params) 

R1     = y(1);
R2     = y(2);
R3     = y(3);
M      = y(4);
I      = y(5);
R1star = y(6);
R2star = y(7);
R3star = y(8);
HC8    = y(9);
HC10   = y(10);
C8     = y(11);
bmax   = 1/80*60;

mu   = 0.001; 
kon  = 0.1*60;
koff = 10*60;
k1C8on    = kon;
k1C8off   = koff;
k1HC8on   = params(1);
k1HC8off  = params(2);
k2HC10on  = kon;
k2HC10off = koff;
k2HC8on   = params(1);
k2HC8off  = params(2);
k3HC8on   = kon;
k3HC8off  = koff;
k3HC10on  = params(1);
k3HC10off = params(2);
kin    = params(3);
katt   = params(4);
alpha2 = params(5);
alpha3 = params(6);


alpha1    = bmax/(1+M/kin);
gamma     = bmax/(1+M/kin)*(R1star/(katt+R1star)+R2star/(katt+R2star));
rho       = alpha1;
delta1    = bmax*(R1star/(katt+R1star)+R3star/(katt+R3star));
delta2    = bmax/(1+R3/kin)*1/(1+R1/kin)*R2star/(katt+R2star);
delta3    = bmax/(1+R2/kin)*(R1star/(katt+R1star)+R3star/(katt+R3star));

dR1dt     = alpha1-k1C8on*C8*R1+k1C8off*R1star-k1HC8on*HC8*R1+k1HC8off*R1star-mu*R1;
dR2dt     = alpha2-k2HC10on*HC10*R2+k2HC10off*R2star-k2HC8on*HC8*R2+k2HC8off*R2star-mu*R2;
dR3dt     = alpha3-k3HC8on*HC8*R3+k3HC8off*R3star-k3HC10on*HC10*R3+k3HC10off*R3star    -mu*R3;
dMdt      = gamma-mu*M;
dIdt      = rho-mu*I-(delta1+delta2+delta3)*I;
dR1stardt = k1C8on*C8*R1-k1C8off*R1star+k1HC8on*HC8*R1-k1HC8off*R1star-mu*R1star;
dR2stardt = k2HC10on*HC10*R2-k2HC10off*R2star+k2HC8on*HC8*R2-k2HC8off*R2star-mu*R2star;
dR3stardt = k3HC8on*HC8*R3-k3HC8off*R3star+k3HC10on*HC10*R3-k3HC10off*R3star    -mu*R3star;
dHC8dt    = delta1*I-k1HC8on*HC8*R1+k1HC8off*R1star-k2HC8on*HC8*R2+k2HC8off*R2star-k3HC8on*HC8*R3+k3HC8off*R3star-mu*HC8;
dHC10dt   = delta2*I-k2HC10on*HC10*R2+k2HC10off*R2star -k3HC10on*HC10*R3+k3HC10off*R3star -mu*HC10;
dC8dt     = delta3*I-k1C8on*C8*R1+k1C8off*R1star-mu*C8;

dydt = [dR1dt;dR2dt;dR3dt;dMdt;dIdt;dR1stardt;dR2stardt;dR3stardt;dHC8dt;dHC10dt;dC8dt];

