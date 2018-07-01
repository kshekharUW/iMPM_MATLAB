clc; clear; close all;
% problem inputs
nel=[2,2]; npe=[2,2];
nx=nel+1;
xmin=[0,0];xmax=[1,1];
Tmax =10.5;
% get phase setup info
[ identifier,Number_of_regions,xcenter,box_inner,box_outer,Hwater,v_initial,rho_initial,mu,K,E,lambda,material,gravity ] = ...
    ipTest();
rhoWater=1000;

% initialize the grid
[dx,Nn,Ne,x,icon,numOC,ntop,nbot,nright,nleft,freeDOFs]=initgrid(nx,xmin,xmax);
lx=dx(1); ly=dx(2);

% assemble matrices
[L,C,DxM,DyM,G,M,DivM,detJ]=assemble(dx,Nn,Ne,icon);
gp=1/(sqrt(3))*[-1,1,1,-1;
                -1,-1,1,1];  % Gauss points
            
% iniitalize the material points
[np,xp,xp0,vp,rp,Fp,rhop,rho0p,sigp,pp,mp,vol_dil,MP_tag]= ...
    initpar(npe,dx,icon,x,Ne,Number_of_regions,xcenter,box_inner,box_outer,v_initial,rho_initial,identifier);

% boundary setup
[pSupportDOF,pFreeDOF,vSupportDOF,vFreeDOF,vBC,pressureBC,qn]=bcSetup(Nn,x,gravity,rhoWater,nbot,ntop,nleft,nright);
ntopright=[];

% time stepping
CFL = 0.05;
dt0 = CFL * min(dx(1),dx(2))/1;
t=0; counter = 0; dt=dt0;

% calculate the natural coordinates
[xpn,nep]=natcoords(xp,dx,xmin,nx-1);
% map material point data to the grid
[mv,v]=maptogridIC(Nn,Ne,icon,xpn,dx,dt,nep,np,mp,vp,sigp,pp,rho0p,gravity,material,freeDOFs,nbot,ntop,nleft,nright,ntopright);
nodePressure= zeros(Nn,1);
b=zeros(2*Nn,1); 
vstar=vBC; 
vnew =zeros(2*Nn,1); 

picCounter=0;
% load ic4x4_t200;
% vstar=vnew;
Re=1;
% dtSwitch=0;
t=t-dt;
while(t<=Tmax)
%     if t>dt*100 && dtSwitch==1;
%         dtSwitch=0;
%         t=0; xp=xp0;
%         clf
%     end
    
    t=t+dt; fprintf('time=%f\n',t);    

    %%% Solve pressure Poisson's equation 
%     prhs = -rhoWater/dt * DivM*vstar + qn;
%     nodePressure(pFreeDOF) = L(pFreeDOF, pFreeDOF)\...
%         (prhs(pFreeDOF)-L(pFreeDOF,pSupportDOF)*pressureBC(pSupportDOF));
%     nodePressure(pSupportDOF) = pressureBC(pSupportDOF);
    
    L(pSupportDOF,pSupportDOF)=1e20;
    prhs = -rhoWater/dt * DivM*vstar + qn;  
    prhs(pSupportDOF)=0;
    nodePressure = L\prhs;
    %%% end
    
    %%% Calculate new velocity correcting for pressure gradient and apply BC        
    for node=1:Nn
        if ismember(node,[nleft,nright])
            gpx=0;
        else
            gpx = (nodePressure(node+1)-nodePressure(node-1))/(2*lx);
        end

        if ismember(node,[ntop,nbot])
            gpy=0;
        else
            gpy = (nodePressure(node+nx(1))-nodePressure(node-nx(1)))/(2*ly);
        end
%         vtx = - dt/rhoWater*gpx;
%         vty = - dt/rhoWater*gpy;
        vnew(node)    = vstar(node)    - dt/rhoWater*gpx;
        vnew(node+Nn) = vstar(node+Nn) - dt/rhoWater*gpy;
    end
    

    vnew(vSupportDOF) = vBC(vSupportDOF);   
    %%% end  
   
    %%% plotting and data storage
    storeInterval=dt;    
    if (mod(t-dt,storeInterval)<1e-8 || storeInterval-mod(t-dt,storeInterval)<1e-8)
        picCounter=picCounter+1;
        plotProj( 1,Nn,xp,vp,x,sigp,xmin,xmax,t-dt,nx,nodePressure,vnew,picCounter,nbot,ntop,nleft,nright,0,0 );
        t;
        clf
    end
    %%% end

    %%% Add viscosity forces using gauss pt integration
    mu=rhoWater/Re;
    accn=zeros(2*Nn,1);  
    tau=zeros(2,2);
    for i=1:Ne
        iv=icon(:,i);
        gv=zeros(2,2);
        
        A = (0.5/ly)*(vnew(iv(1)+Nn) - vnew(iv(2)+Nn) + vnew(iv(3)+Nn) - vnew(iv(4)+Nn));
        B = (0.5/lx)*(vnew(iv(1))    - vnew(iv(2))    + vnew(iv(3))    - vnew(iv(4)));
        for j=1:4 % loop over 4 gauss points
            ss=gp(1,j);tt=gp(2,j);
            gradN = getGrad([ss,tt]);
            
            gv(1,1)=2/lx*gradN(1,:)*vnew(iv);   % velocity gradient 
            gv(1,2)=2/ly*gradN(2,:)*vnew(iv);   % velocity gradient 
            gv(2,2)=2/ly*gradN(2,:)*vnew(iv+Nn);% velocity gradient 
            gv(2,1)=2/lx*gradN(1,:)*vnew(iv+Nn);% velocity gradient 
            
            d      = 0.5*(gv+gv');          % rate of strain at gauss point j
            ddev   = d-1/3*trace(d)*eye(2); % rate of shear at gauss point j             
            
            dhat(1,1)= -2/3*A*ss/(lx/2) + 1/3*B*tt/(ly/2); % enhanced velocity contributions
            dhat(2,2)= -2/3*B*tt/(ly/2) + 1/3*A*ss/(lx/2); % enhanced velocity contributions
            dhat(1,2)= 0; dhat(2,1)=0;       % enhanced velocity contributions
            
            tau      = 2*mu*(ddev); 
            fgp      = tau*gradN*lx*ly/4;
            
            accn(iv)   = accn(iv)   - fgp(1,:)';
            accn(iv+Nn)= accn(iv+Nn)- fgp(2,:)';
        end
    end 
   
    
    %%% 

    %%% Map to points
    [xp,vp,particlePressure]=maptopointsPC(Nn,xmin,nx,np,xp,vp,icon,vnew,v,nodePressure,dx,dt,1);
    %%%% end
    
%     %%% Map to grid
%     [xpn,nep]=natcoords(xp,dx,xmin,nx-1);
%     [mv,vnew]=maptogridIC(Nn,Ne,icon,xpn,dx,dt,nep,np,mp,vp,sigp,pp,rho0p,gravity,material,freeDOFs,nbot,ntop,nleft,nright,ntopright);
%     vnew(vSupportDOF) = vBC(vSupportDOF);
%     %%% end
    
    v=vnew;
    vstar=vnew + accn*dt;
end

