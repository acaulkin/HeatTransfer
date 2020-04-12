%Andrew Caulkins
%Convergence Study Code
%___________________________

clc
clear
close all

%Initial Values:
N_y=48;
N_x=2*N_y;
L=1;
Pe=10;

%Lengths:
dx=(3*L)/N_x;                         
dy=L/N_y;

%Spacings:
x = linspace(dx,1-dx, N_x+1);
y = linspace(dy,1-dy, N_y+1);

%Defining Meshes:
[X,Y] = meshgrid(x,y);
[T,Loc_Nus,Q]=Project2_2(N_y,Pe);
x_L=linspace(0,L,length(Loc_Nus));


figure(2)
plot(x_L,Loc_Nus,'go-','linewidth',2)
xlabel('x/L', 'fontsize', 20);
ylabel('Local Nusslet Number', 'fontsize', 20);
grid on;
set(gca, 'fontsize', 16);



Ny=[12,24,48,96]/2;

for i=1:length(Ny)
  [T,Loc_Nus,Q] =Project2_2(Ny(i),Pe);
  x_L=linspace(0,L,length(Loc_Nus));
  figure(3)
plot(x_L,Loc_Nus,'*-')
xlabel('x/L')
ylabel('Local Nusslet Number');
grid on
h= legend('N_x =  6','N_x =  12','N_x = 24','N_x = 48');
title('All Local Nusselt Numbers')
hold on
end

