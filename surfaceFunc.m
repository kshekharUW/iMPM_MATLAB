function [signF,sGrad,pBoundary,signF2] = surfaceFunc( Nn,Ne,icon,nep,xp,np,x,rp,numOC,dx,mv )
%%% map signed function to grid
signF=ones(Nn,1)*10;
% for ip=1:np;
%     i=nep(ip);
%     iv=icon(:,i);
%     for jj=1:4
%         lx=x(iv(jj),:)-xp(:,ip)';
%         dist=sqrt(lx(1)^2+lx(2)^2) - max(dx)*2;
%         if dist<signF(iv(jj))
%             signF(iv(jj))=dist;   
%         end
%     end
% end

for ip=1:np;
    i=nep(ip);
%     iv=icon(:,i);
    for jj=1:Nn
        lx=x(jj,:)-xp(:,ip)';
%         dist=sqrt(lx(1)^2+lx(2)^2) - max(dx);
        dist=sqrt(lx(1)^2+lx(2)^2) - 4*sqrt(dx(1)*dx(2)/16/pi);
        if dist<signF(jj)
            signF(jj)=dist;   
        end
    end
end

signF2=signF;
index0=find(signF<-1e-14);
signF(index0)=-1;


% for ip=1:np;
%     i=nep(ip);
%     iv=icon(:,i);
%     for jj=1:Nn
%         lx=x(jj,:)-xp(:,ip)';
%         signF(jj)=signF(jj)+sqrt(lx(1)^2+lx(2)^2);        
%     end
% end

% gradient of field
nds=[0,1,1,0;
     0,0,1,1];
 sGrad=zeros(Nn,2);
 
 for i=1:Ne
     iv=icon(:,i);
     for j=1:4
         g=shapeg(nds(:,j),dx);
         sGrad(iv(j),:)= sGrad(iv(j),:) + (g'*signF(iv))';
     end
 end
sGrad(:,1)=sGrad(:,1)./numOC;
sGrad(:,2)=sGrad(:,2)./numOC;

% set gradient to zero if mass =0;
index0=find(mv<=1.e-14);
sGrad(index0,1)=0;
sGrad(index0,2)=0;
% cleanup
index0=find(abs(sGrad(:,1))<=1.e-14);
sGrad(index0,1)=0;
index0=find(abs(sGrad(:,2))<=1.e-14);
sGrad(index0,2)=0;

% Determine pressure boundary
pBoundary=[];
for i=1:Nn
    absNorm=norm(sGrad(i,:));
    if absNorm>1.e-14
        pBoundary(end+1)=i;
    end
end



end

