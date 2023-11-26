%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB code golden_funct1.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alpha1,falpha1] = golden_funct1(x,search)
nmax=200; addx=0.02; adddx_mult=1.2;
epsilon1 = 1.0e-4; epsilon2 = 1.0e-4; % convergence criteria
a = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given initial pt, find the range of sol
fa=func_ALM(x+a*search);n=1;
b=a+addx; fb=func_ALM(x+b*search);
   while (fb < fa),
     fa=fb;
     addx=addx*adddx_mult; b=b+addx;
     fb=func_ALM(x+b*search);
     n=n+1;
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau = (3-sqrt(5))/2;
alpha1 = a*(1-tau) + b*tau;
alpha2 = a*tau + b*(1-tau);
falpha1 = func_ALM(x+alpha1*search);
falpha2 = func_ALM(x+alpha2*search);
for i= 1:nmax
    if falpha1 > falpha2
     a = alpha1;
     alpha1 = alpha2;
     falpha1 = falpha2;
     alpha2 = tau*a + (1-tau)*b;
     falpha2 = func_ALM(x+alpha2*search);
    else
     b = alpha2;
     alpha2 = alpha1;
     falpha2 = falpha1;
     alpha1 = tau*b + (1-tau)*a;
     falpha1 = func_ALM(x+alpha1*search);
    end
 
   if ((abs(falpha1-falpha2)< epsilon1) || (abs(alpha1-alpha2)< epsilon1)),
   break;
   end
end


