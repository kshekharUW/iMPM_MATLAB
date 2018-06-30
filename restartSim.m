function [counter,t,dt,xp,vp,sigp] = restartSim(name,dt,timescaler)
counter=str2num(name(10:end));
t=counter*dt;
dt=dt/timescaler;

M = csvread(strcat('data\',name));

xp=M(:,[1,2])';
vp=M(:,[3,4])';
for i=1:length(xp(1,:))
    sigp(:,i)= [M(i,5);M(i,5);0;0];   
end
