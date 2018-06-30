function [np,xp,xp0,vp,rp,Fp,rhop,rho0p,sigp,pp,mp,vol_dil,MP_tag]=initpar(npe,dx,icon,x,Ne,...
    nrg,x_at_center,radius_inner,radius_outer,vin,rhoin,identifier)
np=0;
wt=1./npe;
pvol=prod(dx./npe);

for i=1:Ne;
    ii=icon(:,i);
    for kr=1:nrg;
        
        r_inner = radius_inner(kr,:);
        r_outer = radius_outer(kr,:);
        xcenter = x_at_center(kr,:);
        for k1=1:npe(1);
            for k2=1:npe(2);
                xi=[k1-0.5,k2-0.5].*wt;
                n=shape(xi);
                xtest=[(x(ii,1)'*n);(x(ii,2)'*n)];
%                 if( (sum((xtest-xcenter(kr,:)').^2)-rout(:,kr)^2) <=0 ...
%                         && (sum((xtest-xcenter(kr,:)').^2)-rin(:,kr)^2) >=0 )

                if(xtest(1) >= xcenter(1)-r_outer(1) && xtest(1) <= xcenter(1)+r_outer(1) ...
                        && xtest(2) >= xcenter(2)-r_outer(2) && xtest(2) <= xcenter(2)+r_outer(2)  )
%                     if xtest(2)+0.8*xtest(1)-2<0
                    if 1
                        np=np+1;
                        MP_tag(np) = identifier(kr);
                        xp(:,np)=xtest;
                        xp0(:,np)=xtest;
                        vp(:,np)= vin(kr,:);
                        rhop(np)=rhoin(kr);
                        rho0p(np)=rhop(np);
    %                     sigp(:,np)=[0;0;0;0];
                        pressure=-(2-xp(2,np))*1000*10;
                        sigp(:,np)=[pressure;pressure;0;0];
                        pp(np)=pressure;
    %                     sigp(:,np)=[0;0;0;0];
                        mp(np)=pvol*rhop(np);
                        Fp(:,:,np)=eye(2);
                        vol_dil(np) = 0;
                        rp(np)=sqrt(mp(np)/rhop(np)/pi);
                    end
                end;
            end;
        end;
    end;
end;