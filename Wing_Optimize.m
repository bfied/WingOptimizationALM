%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB code ALM_DFP.m

%fast and reliable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n_of_var -> number of design variables
% A -> approximation of inverse of the hessian matrix for DFP
% LAMBDA, BETA -> Lagrange Multipliers
% RK -> penalty parameter
clc;
clear all;close all
global n_of_eqc n_of_iqc LAMBDA RK BETA EQCON IQCON FVALUE ...
    nd_area rho p_mat S W Velocityi

addpath('airfoils');
FOIL=importselig('mh60.dat'); %function to import airfoils calculate area
plot(FOIL.upper(:,1),FOIL.upper(:,2));hold on;
plot(FOIL.lower(:,1),FOIL.lower(:,2));
title(FOIL.title);axis([0 1 -0.7 0.7]);
%hold off;

Coordinates= FOIL.normalized;
nd_area= FOIL.area;
rho=1.227; %mass of air kg/m3 STP
p_mat= 40; %mass of Depron style foam kg/m3
charlength= .5 
u= 1.85e-5;
Velocityi= 8; %m/s initial guess at velocity
Re= rho*Velocityi*charlength/u;
%Re=1e6;
alpha= [0:0.08:5];
XFOIL=xfoil(Coordinates,alpha,Re,0);


[val,index]=min(XFOIL.CDp); %finds minimum value of angle of attack

Cl= XFOIL.CL(index);    %minimum Cl at specified AoA
CL= 0.95*Cl;             %2d to finite wing correction 
%CL= 0.95*CLc;


W= 0.25*9.81; %kg
AREAREQ=2*W/(rho*Velocityi^2*CL);
S=AREAREQ;
%calculated fixed wing area for airfoil
%{
design variables
   b    Cr    Ct 
_______________________
[ x(1), x(2), x(3) ]

%}
n_of_var = 4;  %number of design variables
x = [1 0.4 0.08 Velocityi+2];     %starting point
n_of_eqc = 1; %number of equality constraints
n_of_iqc = 8; %number of inequality constraints
epsilon1 = 1e-7; epsilon2 = 1e-7; 
epsilon3 = 1e-7;
delx = 1e-4;
nmax=1000;
A = eye(length(x));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RK=1.0;

if (n_of_eqc ~= 0)
LAMBDA = zeros(1,n_of_eqc); 
else
LAMBDA = zeros(1,1); 
end


if (n_of_iqc ~= 0)
BETA = zeros(1,n_of_iqc); checkconstr = zeros(1,n_of_iqc);
else
BETA = zeros(1,1); checkconstr = zeros(1,1);
end



fprintf('Initial function value =  %7.4f\n ',FVALUE)
fprintf(' No.       x-vector       rk          f(x)        |Cons.| \n')
fprintf('__________________________________________________________\n')

falpha_prev = func_ALM(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iALM = 1:nmax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if iALM==1
    deriv = grad_vec(x,delx,n_of_var);
    search = -deriv;
    [alpha,falpha] = golden_funct1(x,search);
    x = x + alpha*search; falpha = func_ALM(x);   
    else
     deriv = grad_vec(x,delx,n_of_var);
     deltax = (alpha*search); deltag = deriv-deriv_prev;
     term1 = (deltax'*deltax)/(deltax*deltag');
     term2 = (A*deltag'*deltag*A)/(deltag*A*deltag');
     A = A + term1 - term2;

     search = -A*deriv'; search = search';
     [alpha,falpha] = golden_funct1(x,search);
     x = x+alpha*search; falpha = func_ALM(x);
    end
    
     checkconstr1 = max([IQCON;checkconstr]);    
     fprintf('%3d  % 8.3d %8.3d % 8.3d % 8.3e %8.3e %8.3e \n',iALM,x,RK,FVALUE) %norm([EQCON checkconstr1]))
        
      %if abs(falpha-falpha_prev)<epsilon1 || norm(deriv)<epsilon2
      %break;
      %end
        
     falpha_prev = falpha; deriv_prev = deriv;

LAMBDA = LAMBDA + 2*RK*EQCON;
BETA = BETA + 2*RK*(max([IQCON; -BETA./(2*RK)]));

ixdv(iALM,:)=x;ilambda(iALM,:)=LAMBDA; ibeta(iALM,:)=BETA;
if (iALM > 1)
  if  max( abs(ixdv(iALM-1,:)   -ixdv(iALM,:)) ) < epsilon3  && ...
      max( abs(ilambda(iALM-1,:)-ilambda(iALM,:)) ) < epsilon3 && ...
      max( abs(ibeta(iALM-1,:)  -ibeta(iALM,:)) ) < epsilon3 
  break;
  end
end

RK=RK*1.1;                  %edit these values very slowly
if (RK > 2); RK=2; end    %<--- usually not above 2



xx(iALM,:)=x;
fxx(iALM,1)=FVALUE;
gxx(iALM,:)=IQCON;

if (iALM > 1)
 if  max( abs(xx(iALM,:)- xx(iALM,:)) ) < epsilon3 && ... 
        ( abs(fxx (iALM-1,1)-fxx(iALM,1)) < epsilon3) 
 break;
 end
end

end

%{
fprintf('_____________________________________________________\n\n')
fprintf('Lagrange Multipliers(eq & ineq cnst): \n\n')
disp([LAMBDA])
disp([BETA])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}


figure(2)
plot(xx(:,1),'*r-');hold;
plot(xx(:,2),'*b-');
plot(xx(:,3),'or-');
plot(xx(:,4),'ob-');
grid
hold off

figure(3)
plot(fxx(:,1),'*r-');
grid
hold off

figure(4)
plot(gxx(:,1),'ob-'); hold;
plot(gxx(:,2),'ob-');
plot(gxx(:,3),'ob-');
plot(gxx(:,4),'ob-');
plot(gxx(:,5),'ob-');
plot(gxx(:,6),'ob-');plot(gxx(:,7),'ob-');


grid
hold off