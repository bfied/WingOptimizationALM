function mass= masscalc(span,root,tip,nd_area)

global p_mat 

%{
function calculates mass as a function of airfoil geometry
%}

%derived equation for volume for integrating length of wing* nondimensional
%area

% v= integralof nd_area*[(tip-root)*x/span+root) from  (0 to span)
vol= nd_area*span*(root+tip); %cm^3

mass= p_mat*vol; % density g/cm^3 * volume 1-side cm^3



end
