function [ identifier,Number_of_regions,xcenter,box_inner,box_outer,Hwater,v_initial,rho_initial,mu,K,E,lambda,material,gravity ] = ipTest( )
identifier  = [1,0];% Material Identifier :- FLUID = 1  ; SOLID = 0;
Number_of_regions=1;

%Test2 taller
xcenter=[1.0,1.0;
         0.5,0.5];
box_inner=[0.0,0.0;
           0.0,0.0];
box_outer=[1.0,1.0
           0.1,0.1];
Hwater=2.0;
      

v_initial=[0.0,0.0;
          0.0,0.0];
rho_initial=[1000;
             10000];

mu = 0.01;
K = 10000; 
E = 1000000;
lambda = 0.01;
material=[K,mu;
          E,lambda];
gravity = 0*10.0; %m/s^2

end

