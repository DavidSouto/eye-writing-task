%load('hip') 

conds = {'upright','inverted','scrambled'};
t = 1;
for cond = 1:3
    for speedid = 1:3
        subplot(3,3,t)

        l = length(T(subcnt,cond,speedid).dx_lock);

        no_error = T(subcnt,cond,speedid).error == 0;

        h = T(subcnt,cond,speedid).dx_lock_int(:,:);
        v = T(subcnt,cond,speedid).dy_lock_int(:,:);    
        
        title(sprintf('%s, speed %d\n', char(conds(cond)), speedid)); 
        time = (1:length(h))./1000; % in seconds
        if ~isempty(h)
            plot_mean_and_stderr(time,h,'k',0);hold on; plot_mean_and_stderr(time,v,'r',0);
        end
        axis([0 3 -15 15]);

    %    x = (1:length(hipdx)).*(1000/Ctrl.ffreq);

        % plot hipdx
        hold on; plot(TRAJ(cond,speedid,1).t,TRAJ(cond,speedid,1).dx,'-k','LineWidth',1); 
        hold on; plot(TRAJ(cond,speedid,1).t,TRAJ(cond,speedid,1).dy,'-k','LineWidth',1);

        ylabel('velocity [deg/sec]');
        xlabel('time');
        t = t + 1;
    end
end

% hip-traj, interpolate between frames
% htime = [1:length(hipx)*(1000/85)/NINTRP hipx];
% calc_hip_dot = 0;
% if calc_hip_dot == 1,
%     
%     wlk = load('../2D-marker-walker-1000mmS.txt');
%     hip_dot = 10;
% 
%     % scale to the screen
%     [line col] = size(wlk);
% 
%     x_coords_id = 1:2:col-1;
%     y_coords_id = 2:2:col-1;
% 
%     % normalize by x ...
%     max_x_pos = max(max(wlk(:,x_coords_id)));       
% 
%     % -1:1 for x, -0.08 to 0.08 for y
%     wlkn(:,x_coords_id(hip_dot)) = wlk(:,x_coords_id(hip_dot))./max_x_pos; 
%     wlkn(:,y_coords_id(hip_dot)) = -wlk(:,y_coords_id(hip_dot))./max_x_pos;
% 
%     % Scale
%     gfx.Scale = 2000; 
% 
%     % x will span the whole screen
%     hipx = Ctrl.xp2deg.*gfx.Scale.*wlkn(:,x_coords_id(hip_dot));
%     hipy = Ctrl.yp2deg.*gfx.Scale.*wlkn(:,y_coords_id(hip_dot))+Ctrl.xpixel/2;
% 
%     % to remove artifacts, interpolate 10 times
%     hipdx = diff(hipx)*(85);
%     hipdy = diff(hipy)*(85);
% 
%     [a b] = butter(2,60/125,'low');
%     hipdx = filtfilt(a,b,hipdx);
%     hipdy = filtfilt(a,b,hipdy);
%     
%     save('hip','hipdx','hipdy')
% 
% else

%% Sub 1
% % fix
% subplot(211);
% title('pursuit');
% SSUB = 1;
%     
% x = abs(TRACE.POWER(:,1:50)).*2/1000;
%     
% % Average across subs
% hold on; errorbar(f(1:50),mean(x),stderr(x));
% xlabel('Hz');
% ylabel('Normalized Power');
% set(gca,'XTick',0:5:30); 
% axis([0 30 0 (4/1000)]);
% 
% % Plot traces
% x = abs(TRACE.POWER(:,1:50)).*2/1000;
%     
% % Could saccades generate this? Should compare to with saccades and without
% subplot(212)
% 
% % Average across subs
% hold on; plot(TRACE.DY','r'); hold on; plot(mean(TRACE.DY),'b')%,mean(x),stderr(x));
% ylabel('target retinal slip [deg/sec]');
% xlabel('time');
% axis([0 1000 -15 15]);
