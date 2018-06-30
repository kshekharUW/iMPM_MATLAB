function [pSupportDOF,pFreeDOF,vSupportDOF,vFreeDOF,vBC,pressureBC,qn]=bcSetup(Nn,x,gravity,rhoWater,nbot,ntop,nleft,nright)

% %%% Define supported and free DOFs for Hydrostatic
% pSupportDOF = ntop; 
% pFreeDOF = setdiff(1:Nn,ntop);
% 
% vxSuportDOF = [nleft,nright];
% vySupportDOF=[nbot+Nn];
% vimposedDOF = [];
% 
% vSupportDOF = [vxSuportDOF,vySupportDOF,vimposedDOF];
% vFreeDOF = setdiff(1:2*Nn,vSupportDOF);
% 
% 
% % Initialize BCs
% pressureBC=zeros(Nn,1);
% pressureBC(pSupportDOF)=-100;
% vBC = zeros(2*Nn,1);
% vBC(vimposedDOF)=0.1;
% 
% % Calculate Pressure gradient load qn
% qn=zeros(Nn,1);
% lns=length(nbot);
% lbase=x(nbot(end),1)-x(nbot(1),1);
% load=(gravity*lbase*rhoWater)/(2*(lns-1));
% qn(nbot)=2*load;
% qn([nbot(1),nbot(end)])=load;
% %%% end

%%%  cavity Flow pressure on TOP ORIGINAL
% pSupportDOF = ntop(2:end-1);
% pSupportDOF = ntop(1:end);
pSupportDOF = ntop(round(length(ntop)/2));
pFreeDOF = setdiff(1:Nn,pSupportDOF);

vxSuportDOF = [nleft,nright];
vySupportDOF=[ntop+Nn,nbot+Nn];
vimposedDOF = ntop(2:end-1);

vSupportDOF = [vxSuportDOF,vySupportDOF,vimposedDOF];
vFreeDOF = setdiff(1:2*Nn,vSupportDOF);
% end

% Initialize BCs
pressureBC=zeros(Nn,1);
pressureBC(pSupportDOF)=0;
vBC = zeros(2*Nn,1);
vBC(vimposedDOF)=1;
qn=zeros(Nn,1);
%%% end

% %%%  cavity Flow pressure on TOP and BOTTOM velocities
% pSupportDOF = ntop(2:end-1);
% pFreeDOF = setdiff(1:Nn,pSupportDOF);
% 
% vxSuportDOF = [nleft,nright];
% vySupportDOF=[ntop+Nn,nbot+Nn];
% vimposedDOF = [ntop(2:end-1),nbot(2:end-1)];
% 
% vSupportDOF = [vxSuportDOF,vySupportDOF,vimposedDOF];
% vFreeDOF = setdiff(1:2*Nn,vSupportDOF);
% % end
% 
% % Initialize BCs
% pressureBC=zeros(Nn,1);
% pressureBC(pSupportDOF)=0;
% vBC = zeros(2*Nn,1);
% vBC(vimposedDOF)=0.1;
% qn=zeros(Nn,1);
% %%% end

% %%%  cavity Flow pressure on TOP TINKERED
% pSupportDOF = ntop(2:end-1);
% pFreeDOF = setdiff(1:Nn,pSupportDOF);
% 
% vxSuportDOF = [nleft(2:end-1),nright(2:end-1)];
% vySupportDOF=[ntop(2:end-1)+Nn,nbot(2:end-1)+Nn];
% vimposedDOF = ntop;
% 
% vSupportDOF = [vxSuportDOF,vySupportDOF,vimposedDOF];
% vFreeDOF = setdiff(1:2*Nn,vSupportDOF);
% % end
% 
% % Initialize BCs
% pressureBC=zeros(Nn,1);
% pressureBC(pSupportDOF)=0;
% vBC = zeros(2*Nn,1);
% vBC(vimposedDOF)=0.1;
% qn=zeros(Nn,1);
% %%% end

% %  cavity Flow pressure on RIGHT
% pSupportDOF = nright;
% pFreeDOF = setdiff(1:Nn,pSupportDOF);
% 
% vxSuportDOF = [nleft,nright];
% vySupportDOF=[ntop,nbot]+Nn;
% vimposedDOF = ntop(2:end-1);
% 
% vSupportDOF = [vxSuportDOF,vySupportDOF,vimposedDOF];
% vFreeDOF = setdiff(1:2*Nn,vSupportDOF);
% % end
% 
% % Initialize BCs
% pressureBC=zeros(Nn,1);
% pressureBC(pSupportDOF)=0;
% vBC = zeros(2*Nn,1);
% vBC(vimposedDOF)=0.1;
% qn=zeros(Nn,1);
% % end

% %  cavity Flow pressure on TOP velocity on left
% pSupportDOF = ntop;
% pFreeDOF = setdiff(1:Nn,pSupportDOF);
% 
% vxSuportDOF = [nleft,nright];
% vySupportDOF=[ntop,nbot]+Nn;
% vimposedDOF = nleft(2:end-1)+Nn;
% 
% vSupportDOF = [vxSuportDOF,vySupportDOF,vimposedDOF];
% vFreeDOF = setdiff(1:2*Nn,vSupportDOF);
% % end
% 
% % Initialize BCs
% pressureBC=zeros(Nn,1);
% pressureBC(pSupportDOF)=0;
% vBC = zeros(2*Nn,1);
% vBC(vimposedDOF)=0.1;
% qn=zeros(Nn,1);
% % end