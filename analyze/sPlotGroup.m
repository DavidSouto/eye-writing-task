%% Check again the walker oscillation is correct, worrying would be oscillations after it stops ...

include = [1:4 6 9:16 18:23]; % 5 had too few trials 
GT = struct([]);
G = struct([]);

n = 1;
for sub = include
    
    % process matches
    load(['match.' int2str(sub) '.mat'], 'MATCH');

    datafile = [int2str( sub ) '.mat'];
    load(datafile,'S'); % too much data, compression

    % calculate average per sub, save as a_xxx
    for cond = 1:3
        for speed = 1:3
        
            G(cond,speed).m(n) = MATCH(cond,speed).m;

            err_free = S(cond,speed).error == 0;

            GT(cond,speed).gdx_lock(n,:) = ones(1,10000)*NaN;
            GT(cond,speed).gdy_lock(n,:) = ones(1,10000)*NaN;

            l = length(S(cond,speed).dx_lock_int);
            GT(cond,speed).gdx_lock(n,1:l) = nanmean(S(cond,speed).dx_lock_int,1);
            GT(cond,speed).gdy_lock(n,1:l) = nanmean(S(cond,speed).dy_lock_int,1);

            fprintf('number errors %d/%d sub %d, cond %d\n', sum(~err_free),numel(err_free), sub, cond);
        end
    end    
    n = n + 1;
end

%% PLOT
figure;
conds = {'upright','inverted','scrambled'};
t = 1;
for cond = 1:3
    for speedid = 1:3
        subplot(3,3,t)

        exclude = [];
        vex = fExclude(exclude, 19);
        
        h = GT(cond,speedid).gdx_lock(~vex,:);
        v = GT(cond,speedid).gdy_lock(~vex,:);         
        
        title(sprintf('%s, speed %d\n', char(conds(cond)), speedid)); 
        time = (1:length(h))./1000; % in seconds
        plot_mean_and_stderr(time,h,'k',0);hold on; plot_mean_and_stderr(time,v,'r',0);
        axis([0 3 -5 20]);

    %    x = (1:length(hipdx)).*(1000/Ctrl.ffreq);

        % plot hipdx
        hold on; plot(TRAJ(cond,speedid,1).t,TRAJ(cond,speedid).dx,'-k','LineWidth',1); 
        hold on; plot(TRAJ(cond,speedid,1).t,TRAJ(cond,speedid).dy,'-k','LineWidth',1);

        ylabel('velocity [deg/sec]');
        xlabel('time');
        t = t + 1;
    end
end
save 'GT';

%% Plot Matches 
for cond = 1:3
    for speed = 1:3
        M(cond,speed) = nanmean(G(cond,speed).m);        
        STERR(cond,speed) = nansem(G(cond,speed).m);
    end
end

%% FIGURE 1
% TBD: 

figure
cols = 'rgb';
for n = 1:3
    hold on; errorbar(M(n,:),STERR(n,:),cols(n));
%    hold on; plot(G(n,:).m,'or');
end
