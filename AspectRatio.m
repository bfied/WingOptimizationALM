function retval= AspectRatio(span,root,tip)

%span = input is wing span
%root = length of airfoil wing root
%tip  = length of airfoil wing tip
%%Area of trapezoidal wing
topArea= span*(root+tip);

%%
AR=(span^2)/topArea;
retval= AR;
end
