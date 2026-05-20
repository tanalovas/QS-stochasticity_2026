
function dydt = sistemODE(~,y,params) 

phiR1   = y(1);
phiR2   = y(2);
phiR3   = y(3);
phiM    = y(4);
phiI    = y(5);
phiR1s  = y(6);
phiR2s  = y(7);
phiR3s  = y(8);
phiHC8  = y(9);
phiHC10 = y(10);
phiC8   = y(11);

% Basic transcription rate 
beta = params(1);
% Molecule degradation rate
mu = params(2);

% kon primary association constant
k1C8on   = params(3);
k2HC10on = params(3);
k3HC8on  = params(3);

% koff primary dissociation constant
k1C8off   = params(4);
k2HC10off = params(4);
k3HC8off  = params(4);

% kon secondary association constant
k1HC8on  = params(5);
k2HC8on  = params(5);
k3HC10on = params(5);

% koff secondary dissociation constant
k1HC8off  = params(6);
k2HC8off  = params(6);
k3HC10off = params(6);

% Inhibition rate
kin = params(7);

% Activation rate
katt = params(8);

% R2 production rate
alpha2 = params(9);

% R3 production rate
alpha3 = params(10);

% Nonlinear rates
alpha1 = beta / (1+phiM/kin);
rho    = alpha1;
gamma  = beta / (1+phiM/kin)*(phiR1s/(katt+phiR1s)+phiR2s/(katt+phiR2s));
delta1 = beta * (phiR1s/(katt+phiR1s)+phiR3s/(katt+phiR3s));
delta2 = beta / (1+phiR3/kin)*1/(1+phiR1/kin)*phiR2s/(katt+phiR2s);
delta3 = beta / (1+phiR2/kin)*(phiR1s/(katt+phiR1s)+phiR3s/(katt+phiR3s));

% Deterministic equations
dphiR1dt   = alpha1 - mu*phiR1 - k1C8on*phiC8*phiR1 + k1C8off*phiR1s - k1HC8on*phiHC8*phiR1 + k1HC8off*phiR1s;
dphiR2dt   = alpha2 - mu*phiR2 - k2HC10on*phiHC10*phiR2 + k2HC10off*phiR2s - k2HC8on*phiHC8*phiR2 + k2HC8off*phiR2s;
dphiR3dt   = alpha3 - mu*phiR3 - k3HC8on*phiHC8*phiR3 + k3HC8off*phiR3s - k3HC10on*phiHC10*phiR3 + k3HC10off*phiR3s;
dphiMdt    = gamma - mu*phiM;
dphiIdt    = rho - mu*phiI - (delta1+delta2+delta3)*phiI;
dphiR1sdt  = k1C8on*phiC8*phiR1 - k1C8off*phiR1s + k1HC8on*phiHC8*phiR1 - k1HC8off*phiR1s - mu*phiR1s;
dphiR2sdt  = k2HC10on*phiHC10*phiR2 - k2HC10off*phiR2s + k2HC8on*phiHC8*phiR2 - k2HC8off*phiR2s - mu*phiR2s;
dphiR3sdt  = k3HC8on*phiHC8*phiR3 - k3HC8off*phiR3s + k3HC10on*phiHC10*phiR3 - k3HC10off*phiR3s - mu*phiR3s;
dphiHC8dt  = delta1*phiI - k1HC8on*phiHC8*phiR1 + k1HC8off*phiR1s - k2HC8on*phiHC8*phiR2 + k2HC8off*phiR2s - k3HC8on*phiHC8*phiR3 + k3HC8off*phiR3s - mu*phiHC8;
dphiHC10dt = delta2*phiI - k2HC10on*phiHC10*phiR2 + k2HC10off*phiR2s - k3HC10on*phiHC10*phiR3 + k3HC10off*phiR3s - mu*phiHC10;
dphiC8dt   = delta3*phiI - k1C8on*phiC8*phiR1 + k1C8off*phiR1s - mu*phiC8;


dydt = [dphiR1dt;dphiR2dt;dphiR3dt;dphiMdt;dphiIdt;dphiR1sdt;dphiR2sdt;dphiR3sdt;dphiHC8dt;dphiHC10dt;dphiC8dt];

end