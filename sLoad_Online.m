%% Load data from eyelink and .log file

Man.currfile = sprintf('%se%ds%db%dt%d.dat',...
util.dir_name, trl.exp, trl.sub, trl.blk, trl.num);

%[TIME, X, Y, Z] = textread(Man.currfile, '%d %f %f %d'); 
[TIME, X, Y, PUP, Z, d, d, d, d] = textread(Man.currfile, '%f %f %f %f %d %f %f %f %d');

X = (X-gfx.scMid(1))./gfx.h_deg2pix;
Y = -(Y-gfx.scMid(2))./gfx.v_deg2pix;

EVE = find(bitand(Z,4));    % events 

errcodes.SAC = 1;
errcodes.BLINK = 2;
errcodes.MSG = 3;
errcodes.SAC_AMP = 4;
errcodes.GAIN = 5;
errcodes.POSITION_CRIT = 6;

SAC = bitget(Z,errcodes.SAC);
BLINK = bitget(Z,errcodes.BLINK);                

synch_time = EVE(1);
stimOn = EVE(3);
stimOff = EVE(4)+10; % stim erased in the next frame
trialOff = EVE(5);

%% Get parameteres of saccades detected by the Eyelink (not always reliable)
SACCADE = [];
if ~isempty(SAC), 

    % beginning (wrong onset at start)
    SACCADE(:,1) = find(diff([0;SAC]) == 1);
                            
    % end (wrong offset at end)
	SACCADE(:,2) = find(diff([SAC;0]) == -1);   
                            
end

% Two-point central difference derivative (see Bahill et al.)
DX = NaN.*ones(length(X),1);
DY = NaN.*ones(length(Y),1);
for n = 11:(length(X)-11); % by 0.02 s
    DX(n) = (X(n+10)-X(n-10))/0.02; 
    DY(n) = (Y(n+10)-Y(n-10))/0.02; 
end

%ZN = fExpandSacc(bitget(Z,errcodes.sacc), SACCADE, 20);

[nsac dummy] = size(SACCADE);
if ~isempty(SACCADE),
    
    % GET X AMPLITUDE
    SACCADE(:,3) = X(SACCADE(:,2))-X(SACCADE(:,1));
    
    % GET Y AMPLITUDE
    SACCADE(:,4) = Y(SACCADE(:,2))-Y(SACCADE(:,1));
    
    % NORM: 
    SACCADE(:,5) = sqrt(SACCADE(:,3).^2 + SACCADE(:,4).^2);

    % Peak Vel
    for t = 1:nsac,
        SACCADE(t,6) = max(abs(DX(SACCADE(t,1):SACCADE(t,2))));

        % Peak Vel DY 
        SACCADE(t,7) = max(abs(DY(SACCADE(t,1):SACCADE(t,2)))); 
    end
    
    % get duration
    SACCADE(:,8) = SACCADE(:,2)-SACCADE(:,1);
end

SSACCADE = SACCADE;

%% expar. crits
% (assuming 1000 hz)
expar.saccade_crit = 1;

expar.start = stimOn; % start quality check from cohere onset and a little before
expar.end = stimOff;

% CHECK SACC, BLINKS, NO SIGNAL DURING CRITICAL INTERVAL
if sum(BLINK(stimOn:stimOff)), 
    
    fprintf('blink error during critical interval\n');
    trl.err = bitset(trl.err,errcodes.BLINK); 
    
end

if sum((SAC > expar.start) & (SAC < expar.end)), % if there is a saccade within (or straddling) the interval
    
    for n=1:length(SACCADE(:,1)),            % saccades length                        
        fprintf('%.2f degrees sacc (crit %.2f)\n',SACCADE(n,5), expar.saccade_crit);
        if SACCADE(n,5)>expar.saccade_crit,
            trl.err = bitset(trl.err,errcode.SAC_AMP); 
        end
    end
    
end

% Pursuit crit
if sum(Z), 
    
    % no interpolation, traces marked as NaN for the extended duration of the saccade
    % if at least one saccade episode    
    %fExpandSacc(bitget(Z,errcodes.sacc), SACCADE, 20);
    SSAC = fExpandSacc(bitget(Z,errcodes.SAC),SACCADE,20);
    DX(logical(SSAC)) = NaN;            

end
                       
trl.gain = nanmean(DX(stimOn:stimOff)./(trl.DIR(trl.pursuit_dir+1).*gfx.pvel));

if ~trl.fixation,
    
    if abs(trl.gain) < 0.8,
        trl.err = bitset(trl.err, errcodes.GAIN);
        fprintf('gain too low: %.2f\n', abs(trl.gain));
    end
    
    trl.rs = abs(nanmean(DX(stimOn:stimOff))) - gfx.pvel;% pvel: pursuit velocity in degrees/s (7 d/s)

end

% in milliseconds
%TIME = 0:(1000*gfx.ifi):length(trl.traj);

traj = @(range,step,mid,dir) dir.*(0:step:range)-dir.*(fix(range*.5))+mid;
    
% displacement per millisecond (PixPerFrame/MsPerFrame)
trl.traj = traj(gfx.range,gfx.pursuit_vel_pix/(1000*gfx.ifi),gfx.scMid(1),trl.DIR(trl.pursuit_dir+1));
    
% ff = 1;
% for tt=1:gfx.trlframes,
%         if (tt > (gfx.trlframes/2-gfx.frames/2)) && (ff <=gfx.frames), % mid traj
%             
%             ff = ff + 1;
%             
%         end
% end    
midTraj = fix(length(trl.traj)./2);
trl.trajOn = midTraj-100;
trl.trajOff = trl.trajOn+length(stimOn:stimOff)-1;
trl.perr = nanmean(X(stimOn:stimOff) - [(trl.traj(trl.trajOn:trl.trajOff)-gfx.scMid(1))./gfx.h_deg2pix]');

%sca
%keyboard
if sum(abs(trl.perr)>2),

    trl.err = bitset(trl.err, errcodes.POSITION_CRIT);

end

%sca
%keyboard