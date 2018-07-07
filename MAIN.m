clc; clear; close all;
% problem inputs
nel=[8,8]; npe=[2,2];
nx=nel+1;
xmin=[0,0];xmax=[1,1];
Tmax =100.0;
% get phase setup info
[ identifier,Number_of_regions,xcenter,box_inner,box_outer,Hwater,v_initial,rho_initial,mu,K,E,lambda,material,gravity ] = ...
    ipTest();
rhoWater=1000;

% initialize the grid
[dx,Nn,Ne,x,icon,numOC,ntop,nbot,nright,nleft,freeDOFs]=initgrid(nx,xmin,xmax);
lx=dx(1); ly=dx(2);

% assemble matrices
[L,C,DxM,DyM,G,M,DivM,detJ]=assemble(dx,Nn,Ne,icon);
            
% iniitalize the material points

[np,xp,xp0,vp,rp,Fp,rhop,rho0p,sigp,pp,mp,vol_dil,MP_tag]= ...
    initpar(npe,dx,icon,x,Ne,Number_of_regions,xcenter,box_inner,box_outer,v_initial,rho_initial,identifier);

% boundary setup
[pSupportDOF,pFreeDOF,vSupportDOF,vFreeDOF,vBC,pressureBC,qn]=bcSetup(Nn,x,gravity,rhoWater,nbot,ntop,nleft,nright);
ntopright=[];

% time stepping
CFL = 0.1;
dt = CFL * min(lx,ly)/1;
t=0;
Tmax=Tmax+dt;

% Flow parameters
Re=1;
mu=rhoWater/Re;

% calculate the natural coordinates
[xpn,nep]=natcoords(xp,dx,xmin,nx-1);
% map material point data to the grid
[mv,v]=maptogridIC(Nn,Ne,icon,xpn,dx,dt,nep,np,mp,vp,sigp,pp,rho0p,gravity,material,freeDOFs,nbot,ntop,nleft,nright,ntopright);
 
% Fictitious time steps to determine initial state
[nodePressure]=solveP(L,DivM,vBC,qn,pSupportDOF,rhoWater,1.0);
[vnew]=solveVtilda(nodePressure,Nn,nbot,ntop,nleft,nright,nx,lx,ly,rhoWater,1.0);
vnew(vSupportDOF) = vBC(vSupportDOF);
%

picCounter=0;
while(t<=Tmax)    
    fprintf('time=%f\n',t);   
    
    %%% Add viscosity forces using gauss pt integration
    [force]=getForce(mu,Nn,Ne,icon,lx,ly,vnew,mv);
    accn = force./[mv;mv];    
    vstar=vnew + accn*dt;
    vstar(vSupportDOF) = vBC(vSupportDOF);
    %%% 

    %%% Solve pressure Poisson's equation 
    [nodePressure]=solveP(L,DivM,vstar,qn,pSupportDOF,rhoWater,dt);
    %%% end
    
    %%% Calculate new velocity correcting for pressure gradient and apply BC        
    [vtilda]=solveVtilda(nodePressure,Nn,nbot,ntop,nleft,nright,nx,lx,ly,rhoWater,dt);
    vnew = vstar + vtilda;
    vnew(vSupportDOF) = vBC(vSupportDOF);   
    %%% end  
   
    %%% plotting and data storage
    storeInterval=dt;    
    if (mod(t-dt,storeInterval)<1e-8 || storeInterval-mod(t-dt,storeInterval)<1e-8)
        plotProj( 1,Nn,xp,vp,x,sigp,xmin,xmax,t-dt,nx,nodePressure,vnew,picCounter,nbot,ntop,nleft,nright,0,0 );
        picCounter=picCounter+1;
    end
    %%% end

%     %%% Map to points
%     [xp,vp,particlePressure]=maptopointsPC(Nn,xmin,nx,np,xp,vp,icon,vnew,v,nodePressure,dx,dt,1);
%     %%%% end
    
%     %%% Map to grid
%     [xpn,nep]=natcoords(xp,dx,xmin,nx-1);
%     [mv,vnew]=maptogridIC(Nn,Ne,icon,xpn,dx,dt,nep,np,mp,vp,sigp,pp,rho0p,gravity,material,freeDOFs,nbot,ntop,nleft,nright,ntopright);
%     vnew(vSupportDOF) = vBC(vSupportDOF);
%     %%% end
  
    v=vnew;    
    t=t+dt;
end

