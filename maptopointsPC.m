function [xp,vp,particlePressure]=maptopointsPC(Nn,xmin,nx,np,xp,vp,icon,vnew,v,nodePressure,dx,dt, timeInt)
% get natural coordinates of particles within cell
[xpn,nep]=natcoords(xp,dx,xmin,nx-1);

vnew = [vnew(1:Nn,:),vnew(Nn+1:2*Nn,:)];
v = [v(1:Nn,:),v(Nn+1:2*Nn,:)];
particlePressure = zeros(np,1);

ap=zeros(2,np); xpF = zeros(2,np);
k1=zeros(2,np); k2=zeros(2,np); k3=zeros(2,np); k4=zeros(2,np);

if timeInt==1
%%% RK4
for ip=1:np
    i=nep(ip);
    iv=icon(:,i);
    n=shape(xpn(:,ip));
      
    k1(1,ip) = v(iv,1)'*n; 
    k1(2,ip) = v(iv,2)'*n;
    
    xpF(1,ip) = xp(1,ip) +  k1(1,ip) *dt/2;
    xpF(2,ip) = xp(2,ip) +  k1(2,ip) *dt/2;    
    [xpn2,nep2] = natcoords(xpF(:,ip),dx,xmin,nx-1);
    i2 = nep2(1); iv2 = icon(:,i2); nn = shape(xpn2(:));
    k2(1,ip) = 0.5*(vnew(iv2,1)+v(iv2,1))'*nn;
    k2(2,ip) = 0.5*(vnew(iv2,2)+v(iv2,2))'*nn;
    
    
    xpF(1,ip) = xp(1,ip) +  k2(1,ip) *dt/2;
    xpF(2,ip) = xp(2,ip) +  k2(2,ip) *dt/2;    
    [xpn2,nep2] = natcoords(xpF(:,ip),dx,xmin,nx-1);
    i2 = nep2(1); iv2 = icon(:,i2); nn = shape(xpn2(:));    
    k3(1,ip) = 0.5*(vnew(iv2,1)+v(iv2,1))'*nn;
    k3(2,ip) = 0.5*(vnew(iv2,2)+v(iv2,2))'*nn;
    
    
    xpF(1,ip) = xp(1,ip) +  k3(1,ip) *dt;
    xpF(2,ip) = xp(2,ip) +  k3(2,ip) *dt;
    [xpn2,nep2] = natcoords(xpF(:,ip),dx,xmin,nx-1);
    i2 = nep2(1); iv2 = icon(:,i2); nn = shape(xpn2(:));  
    k4(1,ip) = vnew(iv2,1)'*nn;
    k4(2,ip) = vnew(iv2,2)'*nn;
    
    xp(1,ip) = xp(1,ip) + dt*(k1(1,ip)/6 + k2(1,ip)/3 + k3(1,ip)/3 + k4(1,ip)/6);
    xp(2,ip) = xp(2,ip) + dt*(k1(2,ip)/6 + k2(2,ip)/3 + k3(2,ip)/3 + k4(2,ip)/6);
    
    %%% update velocity
%     vp(1,ip)=vp(1,ip)+(vnew(iv,1)-v(iv,1))'*n;
%     vp(2,ip)=vp(2,ip)+(vnew(iv,2)-v(iv,2))'*n;      
    vp(1,ip)=vnew(iv,1)'*n;
    vp(2,ip)=vnew(iv,2)'*n; 
    
    %%% Update stress
    nPressure = nodePressure(iv);
    particlePressure(np) = nPressure'*n;    

end
end


if timeInt==2
%%% RK23
for ip=1:np
    i=nep(ip);
    iv=icon(:,i);
    n=shape(xpn(:,ip));
      
    k1(1,ip) = v(iv,1)'*n; 
    k1(2,ip) = v(iv,2)'*n;
    
    xpF(1,ip) = xp(1,ip) +  k1(1,ip) *dt/2;
    xpF(2,ip) = xp(2,ip) +  k1(2,ip) *dt/2;    
    [xpn2,nep2] = natcoords(xpF(:,ip),dx,xmin,nx-1);
    i2 = nep2(1); iv2 = icon(:,i2); nn = shape(xpn2(:));
    k2(1,ip) = 0.5*(vnew(iv2,1)+v(iv2,1))'*nn;
    k2(2,ip) = 0.5*(vnew(iv2,2)+v(iv2,2))'*nn;
    
    
    xpF(1,ip) = xp(1,ip) +  k2(1,ip) *dt*0.75;
    xpF(2,ip) = xp(2,ip) +  k2(2,ip) *dt*0.75;    
    [xpn2,nep2] = natcoords(xpF(:,ip),dx,xmin,nx-1);
    i2 = nep2(1); iv2 = icon(:,i2); nn = shape(xpn2(:));    
    k3(1,ip) = 0.75*(vnew(iv2,1)+v(iv2,1))'*nn;
    k3(2,ip) = 0.75*(vnew(iv2,2)+v(iv2,2))'*nn;
    
    xp(1,ip) = xp(1,ip) + dt*( 2/9*k1(1,ip) + 1/3*k2(1,ip) + 4/9*k3(1,ip) );
    xp(2,ip) = xp(2,ip) + dt*( 2/9*k1(2,ip) + 1/3*k2(2,ip) + 4/9*k3(2,ip) );
    
    xpF(1,ip) = xp(1,ip);
    xpF(2,ip) = xp(2,ip);
    [xpn2,nep2] = natcoords(xpF(:,ip),dx,xmin,nx-1);
    i2 = nep2(1); iv2 = icon(:,i2); nn = shape(xpn2(:));  
    k4(1,ip) = vnew(iv2,1)'*nn;
    k4(2,ip) = vnew(iv2,2)'*nn;    
end
end


if timeInt==0
%%% Euler
for ip=1:np;
    i=nep(ip);
    iv=icon(:,i);    
    n=shape(xpn(:,ip));
    g=shapeg(xpn(:,ip),dx);  
    
    %%% update velocity
    vp(1,ip)=vp(1,ip)+(vnew(iv,1)-v(iv,1))'*n;
    vp(2,ip)=vp(2,ip)+(vnew(iv,2)-v(iv,2))'*n;
    
    %%% update position
    xp(1,ip)=xp(1,ip)+0.5*dt*(vnew(iv,1)+v(iv,1))'*n;
    xp(2,ip)=xp(2,ip)+0.5*dt*(vnew(iv,2)+v(iv,2))'*n;
   
%     %%% Update stress
%     nPressure = nodePressure(iv);
%     particlePressure(np) = nPressure'*n;
    
%     %%%% Enhanced velocity correction
%     xpnNew= xpn(:,ip)*2-1;
%     Ddivx = 0.25*(vnew(iv(1),2) - vnew(iv(2),2) + vnew(iv(3),2) - vnew(iv(4),2));
%     Ddivy = 0.25*(vnew(iv(1),1) - vnew(iv(2),1) + vnew(iv(3),1) - vnew(iv(4),1));
%     vpHat = [Ddivx*0.5*(1-xpnNew(1)^2); Ddivy*0.5*(1-xpnNew(2)^2)];
%     vp(:,ip)= vp(:,ip)+vpHat;
    
%     xp(1,ip)=xp(1,ip)+dt*vp(1,ip);
%     xp(2,ip)=xp(2,ip)+dt*vp(2,ip);
end
end