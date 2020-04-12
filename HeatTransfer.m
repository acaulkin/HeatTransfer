%Andrew Caulkins
%__________________________


function [T,Loc_Nus,Q]=Project2_2(N_y,Pe)
%The purpose of this Function is to implement a finite difference method in
%order to "solve" the thermal diffusivity problem by populating a sparse
%matrix of temperature value coefficients at specific nodes.


L=1;                               %Plate length
alpha=1*10^-4;                     %Thermal diffusivity 
                            %Peclet number
U=(alpha*Pe)/L;                    %Freestream velocity [m/s]
                            %Size of grid (in y direction)
N_x=2*N_y;                             %Size of grid (in x direction)
T_inf=300;                         %Bulk flow temperature [K]
T_plate=300;                       %Plate temp [K]



dx=(3*L)/N_x;                      %Length between nodes
dy=L/N_y;
x = linspace(dx,1-dx, N_x);        %Position of nodes
y = linspace(dy,1-dy, N_y);

%Defining Sparse Matrix as given in problem statement:
N=(N_y+1)*(N_x+1);
A=sparse(N,N);
b=zeros(N,1);

%Plate Boundary:
space=L/dx;
plate_l= space;
plate_u= 2*space +1;

L = 0;

%Nested for Loop to populate Sparse Matrix:
for j=1:(N_y+1)
for i=1:(N_x+1)
ii  = (j-1)*(N_x+1)+i; %along diagnol of sparse matrix

%Boundary Conditions
%_____________________

%Far Left:
if (i==1)&&(1<=j) && (j<=(N_y+1))           
    A(ii,ii)=1;
    b(ii,1)=T_inf;
end 

if (j==(N_y+1)) && (plate_l<i) && (i<=plate_u)  %Plate
    A(ii,ii)=1;
    b(ii,1)=T_plate; 
    L = L+1;
end
%Top: (analytically solved)
if (j==1) && (1<i) && (i<=(N_x+1))          
    A(ii,ii)=(3/(2*dy));
    A(ii,ii+(N_x+1))=-2/dy;
    A(ii,ii+(2*(N_x+1)))=1/(2*dy);    
end
%Far right Boundary: (analytically solved)
if (i==(N_x+1)) && (1<j) && (j<(N_y+1))    
    A(ii,ii)=(3*U/(2*dx));
    A(ii,ii-1)=-2*U/dx;
    A(ii,ii-2)=U/(2*dx);
end
%Bottom Boundary right of plate: (analytically solved)
if (j==(N_y+1)) && (1<i) && (i<=plate_l)     
    A(ii,ii)=(-3/(2*dy));
    A(ii,ii-(N_x+1))=2/dy;
    A(ii,ii-(2*(N_x+1)))=-1/(2*dy);
end
%Bottom Boundary left of plate: (analytically solved)
if (j==(N_y+1)) && (plate_u<i) && (i<=(N_x+1))     
    A(ii,ii)=(-3/(2*dy));
    A(ii,ii-(N_x+1))=2/dy;
    A(ii,ii-(2*(N_x+1)))=-1/(2*dy);
end
%___________________________________

%Inner Matrix:
if (i>2)&&(j>1)&&(i<(N_x+1))&&(j<(N_y+1))
    A(ii, ii-1) =((-U*2)/dx - alpha/dx^2);
    A(ii , ii+1) =-alpha/dx^2;
    A(ii, ii-2) =  U/(2*dx);
    A(ii, ii-(N_x+1)) =-alpha/dy^2; 
    A(ii, ii + (N_x+1)) = -alpha/dy^2; 
    A(ii,ii)=((U*3)/(2*dx) - alpha*(-2/dx^2  - 2/dy^2));
end

%Special case:
if (i==2)&&(j>1)&&(j<(N_y+1))
A (ii, ii-1) =((-U*2)/dx - alpha/dx^2); 
A(ii, ii+1) =-alpha/dx^2; 
A(ii, ii-(N_x+1)) = -alpha/dy^2; 
A(ii, ii + (N_x+1)) =-alpha/dy^2; 
A(ii,ii)=((U*3)/(2*dx) - alpha*(-2/dx^2  - 2/dy^2));

b(ii, 1)= - (U/(2*dx))*T_inf; 
end

end

end

T = A\b;
s=full(A);

%Plotting Temperature Contours:
Finish = vec2mat(T,N_x+1);
Finish_1 = flipud(Finish);
figure(1)
contourf(Finish_1)
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/100)     
yt = get(gca, 'YTick');                                 % 'XTick' Values
set(gca, 'YTick', yt, 'YTickLabel', yt/100)   
title('Temperature Contours for Pe Number of 100')
ylabel('y')
xlabel('x')
Pe = 10;

%Calculating Heat Flux along the plate:
Heat_Flux = alpha*(3*Finish(end,N_x/3+1:2*N_x/3+1) - 4*Finish(end-1,N_x/3+1:2*N_x/3+1) + Finish(end-2,N_x/3+1:2*N_x/3+1))/(2*dy);
Heat_Flux = Heat_Flux';
%Determining Local Nusselt Number from Heat Flux:
Loc_Nus = (Heat_Flux*L)/(alpha*(T_plate-T_inf));

%Plotting Local Nusselt Number:
X_nus = linspace(0,1,numel(Loc_Nus));
figure(3)
plot(X_nus,Loc_Nus)
grid on
title('Local Nusselt Number for Pe = 100.0')
xlabel('X/L')
ylabel('Nusselt Number')
% set(gca, 'XTick', xt, 'XTickLabel', xt)    

%Calculating GLOBAL Nusselt Number:
P = linspace(0,1,size(Heat_Flux,1))
Q = trapz(P,Heat_Flux)/(alpha*(T_plate-T_inf))
end