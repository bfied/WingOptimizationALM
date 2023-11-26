function foil=importselig(file)

%{
This Matlab function reads Selig format dat file
airfoil coordinates files from airfoiltools.com
 
http://airfoiltools.com/search/index

 The Input is the name file you want to study : (for example)
		file='n2414.dat';
		foil=importselig(file);

data.title= character name of airfoil
data.upper= n x 2 array of xy coordinates of upper airfoil section
data.lower= n x 2 array of xy coordinates of lower airfoil section
data.area = dimensionless area of airfoil cross section
            multiply this area by length of wing  front to back to
            calculate dimensioned area of airfoil cross section

%brett fiedler bjfiedle@asu.edu
%}

T=importdata(file);

%creates data title string
A=cell2mat(T.textdata(1));
B=num2str(T.data(1));
foil.title= strjoin({A,B});


%generates array of airfoil double in T.data struct
a=size(T.data);  %checks size of T.data struct for nx2 if else fixes 
if a(2)==2
    
else
T.textdata(1,:)=[]; %strips first row title information
T.data(1,:)=[];     %strips first row title information
T.textdata=str2double(T.textdata); %changes text data to double
T.data= cat(2,T.textdata,T.data); %concatenates textdata and data columns
end

%splits T.data into useable x,y output arrays for upper and lower

Tsize=size(T.data);

foil.normalized= T.data;
foil.upper= flip(T.data(1:(Tsize(1)-1)/2+1,:)); 
foil.lower= T.data((Tsize(1)-1)/2+1:Tsize,:);

%data.upper= flip(T.data(1:(Tsize(1))/2,:));  test import file
%data.lower= T.data((Tsize(1))/2:Tsize,:);

%calculates upper and lower areas with trapezoidial area function
Area_U=trapz(foil.upper(:,1),foil.upper(:,2));
Area_L= abs(trapz(foil.lower(:,1),foil.lower(:,2)));

foil.area= Area_U+ Area_L;
%foil.xfoil= concatat
end


 
