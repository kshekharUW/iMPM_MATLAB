function [v]=damping(v,icon,Ne,Nn,alpha)
v = [v(1:Nn,:),v(Nn+1:2*Nn,:)];
modes=[1,-1,-1,1;
       1,1,-1,-1;
       1,-1,1,-1;
      -1,1,1,-1];
 
% modes2=[1,-1,-1,1;
%        1,1,-1,-1;
%        1,-1,1,-1;
%       -1,-1,1,1]; 
%   
% modesx=[-1,1,1,-1;
%          1,-1,-1,1;
%          -1,-1,1,1;
%         -1,1,-1,1;
%          0,0,0,0];
%  modesy=[-1,1,1,-1;
%           1,-1,-1,1;
%          -1,-1,1,1;
%          -1,1,-1,1;
%           0,0,0,0];

% modes=[1,-1,1,-1;
%       -1,1,1,-1;
%        1,-1,1,-1;
%       -1,1,1,-1];

modes=[1,-1,1,-1];

if length(v(1,:))==2
    for el = 1:Ne
       iv=icon(:,el);
       v(iv,1) = v(iv,1) - alpha*(modes(1,:)*v(iv,1)*modes(1,:)');
       v(iv,2) = v(iv,2) - alpha*(modes(1,:)*v(iv,2)*modes(1,:)');
    end
else
    for el = 1:Ne
       iv=icon(:,el);
       v(iv,1) = v(iv,1) - alpha*(modes(1,:)*v(iv,1)*modes(1,:)'+...
                                 modes(2,:)*v(iv,1)*modes(2,:)'+...
                                 modes(3,:)*v(iv,1)*modes(3,:)'+...
                                 modes(4,:)*v(iv,1)*modes(4,:)');
    end
end
  
 v=[v(:,1);v(:,2)];
