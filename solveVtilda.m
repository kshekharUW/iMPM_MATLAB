function [vtilda]=solveVtilda(nodePressure,Nn,nbot,ntop,nleft,nright,nx,lx,ly,rhoWater,dt)
vtilda=zeros(2*Nn,1);
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
    vtilda(node)    = - dt/rhoWater*gpx;
    vtilda(node+Nn) = - dt/rhoWater*gpy;
end