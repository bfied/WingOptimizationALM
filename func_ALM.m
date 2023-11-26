function y = func_ALM(x)
global n_of_eqc n_of_iqc LAMBDA RK BETA EQCON IQCON FVALUE ...
    nd_area rho S Velocityi W
%%%%%%%%%%%%%%%%%
%{
design variables
   b    Cr   Ct   Vel
_______________________
[ x(1), x(2), x(3), x(4) ]

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%obj function

G=x(3)/x(2); %Taper Ratio
%MAC= 2/3*x(3)*((1+G+G^2)/(1+G));


mass= masscalc(x(1),x(2),x(3),nd_area);  %function weight call
%W=mass*9.81; % 
%CL= (2*W)/(rho*(15)^2*x(1)*MAC);
%Cl= CL/0.9;

AR= AspectRatio(x(1),x(2),x(3));
%if x(4)<30 %experimentation with x(4) as wing sweep
    OE= 1.78*(1-0.045*AR^0.68)-0.64;
%else
%    OE= 4.61*(1-0.045*AR^0.68)*(cosd(x(4)))^0.15-3.1; %Oswalds Efficiency
%    for wing sweep greater than 30degrees, not applicable here
%end
%y=Cd.CD+(CL^2)/(pi()*OE*(x(1)/MAC));    
y=(W^2)/(0.5*rho*x(4)^2*pi()*x(1)^2*OE); %minimize induced drag
%y=(W^2)/(0.5*rho*Velocity^2*pi()*x(1)^2*OE); %minimize induced drag
h(1)=(x(2)+x(3))*x(1)-S;  %equality constraint on span size
%g(1)= -AR+20;
%g(2)=  AR-40;
g(1)=  x(1)-1.2; %constraint from stress
g(2)=  0.2-x(1);
g(3)=  0.15-x(2);
g(4)=   x(2)-0.4;
g(5)=  -x(3)+0.04;
g(6)=   x(3)-x(2);
g(7)=  -mass*9.81+0.18; %g2
g(8)= G-1;

if (n_of_eqc == 0); h(1)=0; end
if (n_of_iqc == 0); g(1)=0; end


EQCON = h; IQCON = g; FVALUE = y;
% Build augmented lagrange function
y = y + ...
    LAMBDA.*EQCON + RK.*EQCON^2 + ...
    sum(BETA.* max([IQCON; -BETA./(2*RK)])) + ...
    sum(RK   *(max([IQCON; -BETA./(2*RK)])).^2);

