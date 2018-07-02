function [force]=getForce(mu,Nn,Ne,icon,lx,ly,vnew,mv)
    gp=1/(sqrt(3))*[-1,1,1,-1;
                    -1,-1,1,1];  % Gauss points

    force=zeros(2*Nn,1);
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
            
%             accn(iv)   = accn(iv)   - fgp(1,:)'./mv(iv);
%             accn(iv+Nn)= accn(iv+Nn)- fgp(2,:)'./mv(iv);
            force(iv)   = force(iv)   - fgp(1,:)';
            force(iv+Nn)= force(iv+Nn)- fgp(2,:)';            
        end
    end