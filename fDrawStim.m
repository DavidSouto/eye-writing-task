% Get response here
function trl = fDrawStim(win,gfx,trl)

     if trl.Eyelink,   
         
            Eyelink('command','begin_realtime_mode(100)');
            Eyelink('Message','TRIAL ON');   % EVENT 1
            
     end
       
    % Define the trajectory
    traj = @(range,step,mid,dir) dir.*(0:step:range)-dir.*(fix(range*.5))+mid;
    
    trl.traj = traj(gfx.range,gfx.pursuit_vel_pix,gfx.scMid(1),trl.DIR(trl.pursuit_dir+1));
%    trl.traj = traj(gfx.range,gfx.pursuit_vel_pix/(1000*gfx.ifi),gfx.scMid(1),trl.DIR(trl.pursuit_dir+1));

    gfx.trlframes = length(trl.traj);
    
    % Start trial
    Screen('Flip', win);        

    gfx.speed_pix_frame = (trl.vel./gfx.h_pix2deg)./gfx.monitor_freq;
%    fprintf('Velocity of cloud %2.2f (%2.2f deg), of target %2.2f (%2.2f deg)\n', gfx.speed_pix_frame, gfx.speed_pix_frame*gfx.h_pix2deg*gfx.monitor_freq, gfx.pursuit_vel_pix, gfx.pursuit_vel_pix*gfx.h_pix2deg*gfx.monitor_freq);
    
    gfx = DotMovement(gfx);

    gfx.sf.x1 = gfx.sf.x;
    gfx.sf.y1 = gfx.sf.y;
    
    if ~trl.stair,
        
        gfx.speed_pix_frame = (-trl.vel./gfx.h_pix2deg)./gfx.monitor_freq;
        %    fprintf('Velocity of cloud %2.2f (%2.2f deg), of target %2.2f (%2.2f deg)\n', gfx.speed_pix_frame, gfx.speed_pix_frame*gfx.h_pix2deg*gfx.monitor_freq, gfx.pursuit_vel_pix, gfx.pursuit_vel_pix*gfx.h_pix2deg*gfx.monitor_freq);

        gfx = DotMovement(gfx);
        gfx.sf.x2 = gfx.sf.x;
        gfx.sf.y2 = gfx.sf.y;

    end
    
    % back to black  for 500 ms
    Screen('DrawDots', win, [0 0], gfx.fixsize, [255 255 255], [trl.traj(1) gfx.scMid(2)], 1);
    Screen('Flip', win); 
    WaitSecs(.5);
    
    % 200 ms, centered on screen centre
    ff = 1;
    for tt=1:gfx.trlframes,

        if (tt >= (gfx.trlframes/2-gfx.frames/2)) && (tt < (gfx.trlframes/2+gfx.frames/2)), % mid traj
            
            % MVT ONSET
            if (trl.Eyelink && ff==1), Eyelink('Message','MOVEMENT ON'); end % event 2

            if ~trl.stair, % in the calib condition display motion in the opposite direction in the lvf or uvf
                                 % the control remains the same as exp1

                    Screen('DrawDots', win, cat(2,gfx.sf.x1(:,ff),gfx.sf.y1(:,ff))', gfx.dotSize, kron([255 255 255],gfx.c(:,ff))', ...
                                [gfx.scMid(1)-gfx.width_pix*.5 (gfx.scMid(2)-gfx.height_pix*.5-trl.ecc)], 1);

                    Screen('DrawDots', win, cat(2,gfx.sf.x2(:,ff),gfx.sf.y2(:,ff))', gfx.dotSize, kron([255 255 255],gfx.c(:,ff))', ...
                                [gfx.scMid(1)-gfx.width_pix*.5 (gfx.scMid(2)-gfx.height_pix*.5--trl.ecc)], 1);                
       
            else

                    Screen('DrawDots', win, cat(2,gfx.sf.x1(:,ff),gfx.sf.y1(:,ff))', gfx.dotSize, kron([255 255 255],gfx.c(:,ff))', ...
                                [gfx.scMid(1)-gfx.width_pix*.5 (gfx.scMid(2)-gfx.height_pix*.5-trl.ecc)], 1);                

            end
   
            % MVT gfx.offset
            if (trl.Eyelink && ff == gfx.frames), Eyelink('Message','MOVEMENT OFF'); end % event 3
            
            ff = ff + 1;
            
        end

        % gfx.offset is half the trajectory in trlframes
        Screen('DrawDots', win, [0 0], gfx.fixsize, [255 255 255], [trl.traj(tt) gfx.scMid(2)], 1);

        Screen('Flip', win);

    end    
    
    Screen('Flip', win);
    
    if trl.Eyelink      
        
        Eyelink('Message','TRIAL OFF'); % event 4  
        
    end
    
end
function gfx = DotMovement(gfx)
    
    gfx.radius = 10 * gfx.h_deg2pix;
    
    dn = gfx.sf.num;
        
    % Pre-allocate for speed
    % gfx.sf.x = sqrt(rand(dn,1)).*gfx.radius_pix;

    gfx.sf.x = rand(dn,1).* gfx.width_pix;    
    gfx.sf.y = rand(dn,1).* gfx.height_pix;

    gfx.c(:,1) = gaussian_border(gfx.sf.x, gfx.width_pix./2, gfx.radius, gfx.sigma); % sigma is half the square
    
    %[gfx.sf.x(:,1),gfx.sf.y(:,1)]=pol2cart(gfx.sf.a, gfx.sf.r);
    gfx.sf.lifetime = round(rand(dn,1).*gfx.lifetime_frames);
        
    gfx.sf.signal = zeros(dn,1);
    
    % define as signal n dots
    gfx.sf.signal(1:gfx.sf.ncohere,1) = 1;
    
    gfx.sf.xv = gfx.speed_pix_frame; %gfx.speed_pix_frame
    gfx.sf.yv = 0;
    
    for f=2:gfx.frames,
        
        % move dots
        gfx.sf.x(:,f) = gfx.sf.x(:,f-1)+gfx.sf.xv;
        gfx.sf.y(:,f) = gfx.sf.y(:,f-1)-gfx.sf.yv; %has to be negative, otherwise wrong direction
        
        % reposition dead dots
        gfx.sf.lifetime = gfx.sf.lifetime-1;
        dead = gfx.sf.lifetime<1;
        
        % dead ones are replaced at a random location
        gfx.sf.x(dead,f) = rand(sum(dead),1).* gfx.width_pix;    
        gfx.sf.y(dead,f) = rand(sum(dead),1).* gfx.height_pix;

        % new life
        gfx.sf.lifetime(dead) = gfx.lifetime_frames;
        
        % 100% coherent
        % update the direction of the noise/dead ones
        % updateDirection = intersect(noise,dead);
    
        % gfx.sf.d(updateDirection) = rand(length(updateDirection),1).*2.*pi;
        % [gfx.sf.xv(updateDirection), gfx.sf.yv(updateDirection)] = pol2cart(gfx.sf.d(updateDirection), gfx.speed_pix_frame);

        % wrap points
        out = gfx.sf.x(:,f) > gfx.width_pix;
        gfx.sf.x(out,f) = 0;
        
        out = gfx.sf.x(:,f) < 0;
        gfx.sf.x(out,f) = gfx.width_pix;

        out = gfx.sf.y(:,f) > gfx.height_pix;
        gfx.sf.y(out,f) = 0;

        out = gfx.sf.y(:,f) < 0;
        gfx.sf.y(out,f) = gfx.width_pix;

        % horizontal position determines dot contrast
        gfx.c(:,f) = gaussian_border(gfx.sf.x(:,f), gfx.width_pix./2, gfx.radius, gfx.sigma); % sigma is half the square
        
    end
    
end
function y = gaussian_border(x,x0,radius,sigma)

    y = x*0+1; % faster init
   
   % Produce a step with soft edges
   low_edg = x < (x0-radius);
   high_edg = x > (x0+radius);
   
   y(low_edg) = exp(-0.5*((x(low_edg)-(x0-radius)).^2)./(2*sigma.^2));
   y(high_edg) = exp(-0.5*((x(high_edg)-(x0+radius)).^2)./(2*sigma.^2));

end