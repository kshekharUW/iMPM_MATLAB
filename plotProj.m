function [ ] = plotProj( fignum,Nn,xp,vp,x,sigp,xmin,xmax,t,nx,nodePressure,v,counter,nbot,ntop,nleft,nright,savepic,savedata ) 
    v = [v(1:Nn,:),v(Nn+1:2*Nn,:)];
    fig=figure(fignum);  
%     clf(fig)  
      
%     [X,Y]=meshgrid(linspace(0,xmax(1),nx(1)),linspace(0,xmax(2),nx(2)), nx(2));
%     plot(x(:,1),x(:,2),'k+');  hold on;
%     surf(X,Y,reshape(nodePressure,[nx(1),nx(2)])'); hold on;    
%     shading interp; 
%     c=colorbar; c.Label.String = 'Pressure(Pa)'; caxis([-100,100])


    [X,Y]=meshgrid(linspace(0,xmax(1),nx(1)),linspace(0,xmax(2),nx(2)), nx(2));
%     xstart=xmax(1)*rand(nx(2),1);
%     ystart=xmax(2)*rand(nx(1),1);
%     xstart=linspace(xmin(1),xmax(1),nx(2));
%     ystart=linspace(xmin(2),xmax(2),nx(1));
%     xstart=X(1:3:end,1:3:end);
%     ystart=Y(1:3:end,1:3:end);
%     streamline(X,Y,reshape(v(:,1),[nx(1),nx(2)])', reshape(v(:,2),[nx(1),nx(2)])',xstart, ystart); hold on; 


    scale_factor = range(x(:,1))/range(x(:,2))/2;
    quiver(x(:,1),x(:,2),v(:,1)*scale_factor,v(:,2)*scale_factor,'b','linewidth',1,'AutoScale','on'); hold on;
    quiver(xp(1,:),xp(2,:),vp(1,:),vp(2,:),'k','AutoScale','on'); hold on;

%     [C,h]=contourf(reshape(x(:,1),nx)',reshape(x(:,2),nx)',reshape(v(:,1),nx)');  hold on
%     set(h,'LineColor','none');
%     c=colorbar; c.Label.String = 'v_x(m/s)'; caxis([-0.1,0.1]);  
    
%     [C,h]=contourf(reshape(x(:,1),nx)',reshape(x(:,2),nx)',reshape(v(:,2),nx)');  hold on
%     set(h,'LineColor','none');
%     c=colorbar; c.Label.String = 'v_y(m/s)'; caxis([-0.1,0.1]);   

%     [C,h]=contourf(reshape(x(:,1),nx)',reshape(x(:,2),nx)',reshape(nodePressure,nx)');  hold on
%     set(h,'LineColor','none');
%     c=colorbar; c.Label.String = 'pressure (Pa)'; caxis([-2,2]);  
    

    plot(xp(1,:),xp(2,:),'k.','Markersize',1,'linewidth',1); hold on
    plot(x(:,1),x(:,2),'bo','Markersize',5);  hold on
    xlabel('x'); ylabel('y');
    title(strcat('t = ',sprintf('%0.5f',t),' s'));
    xdomain=xmax(1)-xmin(1); ydomain=xmax(2)-xmin(2);    
    xlim([xmin(1)-xdomain/25,xmax(1)+xdomain/25]);
    ylim([xmin(2)-ydomain/25,xmax(2)+ydomain/25]);
    ax=gca;
    ax.XTick=([xmin(1):0.5:xmax(1)]);
    ax.YTick=([xmin(2):0.5:xmax(2)]);
    axis equal;  

    set(gca,'fontsize', 30);
    drawnow;
    if savepic==1
        print(fig,strcat('images\img',num2str(counter)),'-dpng','-r100')
%         print(fig,strcat('images\img',num2str(counter)),'-dpdf')
    end


    % write particle data to file
    if savedata==1
        fname=strcat('data\data.csv.',num2str(counter));
        csvwrite(fname,[xp;vp;sigp(1,:)]');   
        fname=strcat('data\nodeData.csv.',num2str(counter));
        csvwrite(fname,[x';v';nodePressure']');
    end
    hold off;
end

