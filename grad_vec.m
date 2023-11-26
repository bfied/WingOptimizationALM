%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB code grad_vec.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deriv = grad_vec(x,delx,n_of_var)
ideriv=0; % 0 for central difference, 1 or others for forward difference

for i = 1:n_of_var
xvec1 = x; xvec2 = x;
xvec2(i) = x(i) + delx; % change i-th DV only.
   if (ideriv == 0),
    xvec1(i) = x(i) - delx;
    deriv(i) = (func_ALM(xvec2) - func_ALM(xvec1) )/(2*delx);
   else
    deriv(i) = (func_ALM(xvec2) - func_ALM(x) )/(delx);
   end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

