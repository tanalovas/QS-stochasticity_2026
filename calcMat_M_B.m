
function [M,B] = calcMat_M_B(params,Eq)


% Basic transcription rate 
beta = params(1);
% Molecule degradation rate
mu   = params(2);

% Primary association constant (kon)
k1C8on   = params(3);
k2HC10on = params(3);
k3HC8on  = params(3);

% Primary dissociation constant (koff)
k1C8off   = params(4);
k2HC10off = params(4);
k3HC8off  = params(4);

% Secondary association constant (kon)
k1HC8on  = params(5);
k2HC8on  = params(5);
k3HC10on = params(5);

% Secondary dissociation constant (koff)
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


% Concentrations at the equilibrium points (Eq)
phiR1   = Eq(1);
phiR2   = Eq(2);
phiR3   = Eq(3);
phiM    = Eq(4);
phiI    = Eq(5);
phiR1s  = Eq(6);
phiR2s  = Eq(7);
phiR3s  = Eq(8);
phiHC8  = Eq(9);
phiHC10 = Eq(10);
phiC8   = Eq(11);


% Nonlinear rates
alpha1 = beta / (1+phiM/kin);
gamma  = beta / (1+phiM/kin)*(phiR1s/(katt+phiR1s)+phiR2s/(katt+phiR2s));
rho    = alpha1;
delta1 = beta * (phiR1s/(katt+phiR1s)+phiR3s/(katt+phiR3s));
delta2 = beta / (1+phiR3/kin)*1/(1+phiR1/kin)*phiR2s/(katt+phiR2s);
delta3 = beta / (1+phiR2/kin)*(phiR1s/(katt+phiR1s)+phiR3s/(katt+phiR3s));


% Initialize matrices
M = zeros(11);
B = zeros(11);



% --- Non-zero entries of matrix M ---

% R1
M(1,1)  = - mu - k1C8on * phiC8 - k1HC8on * phiHC8;
M(1,4)  = - beta * kin/(kin + phiM)^2;
M(1,6)  = k1C8off + k1HC8off;
M(1,9)  = - k1HC8on * phiR1;
M(1,11) = - k1C8on * phiR1;

% R2
M(2,2)  = - mu - k2HC10on * phiHC10 - k2HC8on * phiHC8;
M(2,7)  = k2HC10off + k2HC8off;
M(2,9)  = - k2HC8on * phiR2;
M(2,10) = - k2HC10on * phiR2;

% R3
M(3,3)  = - mu - k3HC8on * phiHC8 - k3HC10on * phiHC10;
M(3,8)  = k3HC8off + k3HC10off;
M(3,9)  = - k3HC8on * phiR3;
M(3,10) = - k3HC10on * phiR3;

% M 
M(4,4) = -beta * kin / (kin+phiM)^2 * (phiR1s/(katt+phiR1s)+phiR2s/(katt+phiR2s)) - mu;
M(4,6) = katt / (katt + phiR1s)^2 * beta * kin / (kin+phiM);
M(4,7) = katt / (katt + phiR2s)^2 * beta * kin / (kin+phiM);

% I 
M(5,1) = beta * phiI * kin/(kin+phiR3) * kin/(kin+phiR1)^2 * phiR2s/(katt + phiR2s);
M(5,2) = beta * phiI * kin/(kin+phiR2)^2 * (phiR1s/(katt+phiR1s) + phiR3s/(katt + phiR3s));
M(5,3) = beta * phiI * kin/(kin+phiR3)^2 * kin/(kin+phiR1) * phiR2s/(katt + phiR2s);
M(5,4) = -beta * kin/(kin+phiM)^2;
M(5,5) = -mu - (delta1+delta2+delta3);
M(5,6) = -beta * phiI * katt/(katt+phiR1s)^2 - beta * phiI * kin/(kin+phiR2) * katt/(katt+phiR1s)^2;
M(5,7) = -beta * phiI * kin/(kin+phiR3) * kin/(kin+phiR1) * katt/(katt+phiR2s)^2;
M(5,8) = -beta * phiI * katt/(katt+phiR3s)^2 -beta * phiI * kin/(kin+phiR2) * katt/(katt+phiR3s)^2;

% R1*
M(6,1)  = k1C8on * phiC8 + k1HC8on * phiHC8;
M(6,6)  = -mu - k1C8off - k1HC8off;
M(6,9)  = k1HC8on * phiR1;
M(6,11) = k1C8on * phiR1;

% R2*
M(7,2)  = k2HC10on * phiHC10 + k2HC8on * phiHC8;
M(7,7)  = - mu - k2HC10off - k2HC8off;
M(7,9)  = k2HC8on * phiR2;
M(7,10) = k2HC10on * phiR2;

% R3*
M(8,3)  = k3HC8on * phiHC8 + k3HC10on * phiHC10;
M(8,8)  = - mu - k3HC8off - k3HC10off;
M(8,9)  = k3HC8on * phiR3;
M(8,10) = k3HC10on * phiR3;

% HC8
M(9,1) = - k1HC8on * phiHC8;
M(9,2) = - k2HC8on * phiHC8;
M(9,3) = - k3HC8on * phiHC8;
M(9,5) = delta1;
M(9,6) = k1HC8off + beta * phiI * katt/(katt+phiR1s)^2;
M(9,7) = k2HC8off;
M(9,8) = k3HC8off + beta * phiI * katt/(katt+phiR3s)^2;
M(9,9) = - mu - k1HC8on * phiR1 - k2HC8on * phiR2 - k3HC8on * phiR3;

% HC10
M(10,1)  = - beta * phiI * kin/(kin+phiR3) * kin/(kin+phiR1)^2 * phiR2s/(katt+phiR2s);
M(10,2)  = - k2HC10on * phiHC10;
M(10,3)  = - k3HC10on * phiHC10 - beta * phiI * kin/(kin+phiR3)^2 * kin/(kin+phiR1) * phiR2s/(katt+phiR2s);
M(10,5)  = delta2;
M(10,7)  = k2HC10off + beta * phiI * kin/(kin+phiR3) * kin/(kin+phiR1) * katt/(katt+phiR2s)^2;
M(10,8)  = k3HC10off;
M(10,10) = - mu - k2HC10on * phiR2 - k3HC10on * phiR3;

%C8
M(11,1)  = - k1C8on * phiC8;
M(11,2)  = - beta * phiI * kin/(kin+phiR2)^2 * (phiR1s/(katt+phiR1s) + phiR3s/(katt+phiR3s));
M(11,5)  =  delta3;
M(11,6)  = k1C8off + beta * phiI * kin/(kin+phiR2) * katt/(katt+phiR1s)^2;
M(11,8)  = beta * phiI * kin/(kin+phiR2) * katt/(katt+phiR3s)^2;
M(11,11) = - mu - k1C8on * phiR1;



% --- Non-zero entries of matrix B ---

% R1
B(1,1)  = alpha1 + mu * phiR1 + k1C8on * phiC8 * phiR1 +  k1C8off * phiR1s +  k1HC8on * phiHC8 * phiR1 +  k1HC8off * phiR1s;
B(1,6)  =  - k1C8on * phiC8 * phiR1 -  k1C8off * phiR1s -  k1HC8on * phiHC8 * phiR1 -  k1HC8off * phiR1s;
B(1,9)  =  k1HC8on * phiHC8 * phiR1 +  k1HC8off * phiR1s;
B(1,11) = k1C8on * phiC8 * phiR1 +  k1C8off * phiR1s;

% R2
B(2,2)  = alpha2 + mu * phiR2 + k2HC10on * phiHC10 * phiR2 +  k2HC10off * phiR2s +  k2HC8on * phiHC8 * phiR2 +  k2HC8off * phiR2s;
B(2,7)  =  - k2HC10on * phiHC10 * phiR2 -  k2HC10off * phiR2s -  k2HC8on * phiHC8 * phiR2 -  k2HC8off * phiR2s;
B(2,9)  =  k2HC8on * phiHC8 * phiR2 +  k2HC8off * phiR2s;
B(2,10) = k2HC10on * phiHC10 * phiR2 +  k2HC10off * phiR2s;

% R3
B(3,3)  = alpha3 + mu * phiR3 + k3HC8on * phiHC8 * phiR3 +  k3HC8off * phiR3s +  k3HC10on * phiHC10 * phiR3 +  k3HC10off * phiR3s;
B(3,8)  =  - k3HC8on * phiHC8 * phiR3 -  k3HC8off * phiR3s -  k3HC10on * phiHC10 * phiR3 -  k3HC10off * phiR3s;
B(3,9)  =  k3HC8on * phiHC8 * phiR3 +  k3HC8off * phiR3s;
B(3,10) = k3HC10on * phiHC10 * phiR3 +  k3HC10off * phiR3s;

% M
B(4,4) =  gamma + mu * phiM;

% I
B(5,5)  =  rho + mu * phiI + (delta1+delta2+delta3) * phiI;
B(5,9)  =  - delta1 * phiI;
B(5,10) =  - delta2 * phiI;
B(5,11) =  - delta3 * phiI;

% R1*
B(6,6)  =  mu * phiR1s +  k1C8on * phiC8 * phiR1 +  k1C8off * phiR1s +  k1HC8on * phiHC8 * phiR1 +  k1HC8off * phiR1s;
B(6,9)  =  - k1HC8on * phiHC8 * phiR1 -  k1HC8off * phiR1s;
B(6,11) = - k1C8on * phiC8 * phiR1 -  k1C8off * phiR1s;

% R2*
B(7,7)  =  mu * phiR2s +  k2HC10on * phiHC10 * phiR2 +  k2HC10off * phiR2s +  k2HC8on * phiHC8 * phiR2 +  k2HC8off * phiR2s;
B(7,9)  =  - k2HC8on * phiHC8 * phiR2 -  k2HC8off * phiR2s;
B(7,10) = - k2HC10on * phiHC10 * phiR2 -  k2HC10off * phiR2s;

% R3*
B(8,8)  =  mu * phiR3s +  k3HC8on * phiHC8 * phiR3 +  k3HC8off * phiR3s +  k3HC10on * phiHC10 * phiR3 +  k3HC10off * phiR3s;
B(8,9)  =  - k3HC8on * phiHC8 * phiR3 -  k3HC8off * phiR3s;
B(8,10) = - k3HC10on * phiHC10 * phiR3 -  k3HC10off * phiR3s;

% HC8
B(9,9) = mu * phiHC8 +  k1HC8on * phiHC8 * phiR1 +  k1HC8off * phiR1s +  k2HC8on * phiHC8 * phiR2 +  k2HC8off * phiR2s +  k3HC8on * phiHC8 * phiR3 +  k3HC8off * phiR3s + delta1 * phiI;

% HC10
B(10,10) = mu * phiHC10 +  k2HC10on * phiHC10 * phiR2 +  k2HC10off * phiR2s +  k3HC10on * phiHC10 * phiR3 +  k3HC10off * phiR3s +  delta2 * phiI;

% C8
B(11,11) = mu * phiC8 +  k1C8on * phiC8 * phiR1 +  k1C8off * phiR1s +  delta3 * phiI;


% Symmetrize B
for i = 1:11
    for j = 1:i
        B(i,j) = B(j,i);
    end
end


end