function [mv,v]=maptogridIC(Nn,Ne,icon,xpn,dx,dt,nep,np,mp,vp,sigp,pp,rhop,gravity,material,freeDOFs,nbot,ntop,nleft,nright,ntopright)
mu = material(1,2);
% zero arrays
mv=zeros(Nn,1);
pv=zeros(Nn,2);
f=zeros(Nn,2);
fstar=zeros(Nn,2);
b=zeros(Nn,2);
v=zeros(Nn,2);
vstar=zeros(Nn,2);
nodePressure=zeros(Nn,1);
%%% map velocities to grid
for ip=1:np;
    i=nep(ip);
    iv=icon(:,i);
    n=shape(xpn(:,ip));
%     volp=mp(ip)/rhop(ip);
    mv(iv)=mv(iv)+mp(ip)*n;
%     nodePressure(iv)=nodePressure(iv)+mp(ip)*sigp(1,ip)*n;
    pv(iv,:)=pv(iv,:)+mp(ip)*[vp(1,ip)*n,vp(2,ip)*n];
end
index1=find(mv>=1.e-14);
index0=find(mv< 1.e-14);
v(index1,1)=pv(index1,1)./mv(index1);
v(index1,2)=pv(index1,2)./mv(index1);
v(index0,1)=0;
v(index0,2)=0;
%%%
v = [v(:,1);v(:,2)];
