% tbd: use only fExpandEvent

%% Load data from eyelink and .log file

% check limits
if Trl.num > length(DATA(:,Col.trl))
    Trl.num = length(DATA(:,Col.trl));
    msgbox('End of File');
end
if Trl.num < 1
    Trl.num = 1;
    msgbox('Start of File');
end

Man.currfile = sprintf('%se%ds%d/%db%dt%d.txt',...
Man.datpath,Trl.exp,Trl.sub, DATA(Trl.num,Col.sess), DATA(Trl.num,Col.blk),DATA(Trl.num,Col.trl));

%%    
if exist(Man.currfile,'file')
%     TMP = textread(Man.currfile);
%     %[TIME, d, d, d, X, Y, Z] 
%     TIME = TMP(:,1);
%     X = TMP(:,7);
%     Y = TMP(:,8);
%     Z = TMP(:,9);
%     
%     
    TMP = textread(Man.currfile);
    %[TIME, d, d, d, X, Y, Z] 
    TIME = TMP(:,1);
    X = TMP(:,6);
    Y = TMP(:,7);
    E = TMP(:,5); % message events
    Z = TMP(:,9); % saccade events
    Z = bitor(Z,E);
else
    fprintf('error: file does not exist\n');
    fprintf('%s\n', Man.currfile);
    error(Trl.num) = bitset(error(Trl.num),Err.num_events);

    TIME = NaN.*ones(4000,1);
    X = NaN.*ones(4000,1);
    Y = NaN.*ones(4000,1);
    Z = [4;zeros(4000-3,1);4;4];
    
end
X = (X-Ctrl.xpixel/2)*Ctrl.xp2deg;
Y = -(Y-Ctrl.ypixel/2)*Ctrl.yp2deg;

EVE = find(bitand(Z,ev.msg));    % events 

%% EVENTS
synch_time = EVE(2);
fixOn = EVE(3);
stimOn = EVE(4);

if numel(EVE)<5
    stimOff = numel(X);
else
    stimOff = EVE(5);
end
%recOff = EVE(4); % add this 
recOff = stimOff;
    
%% Get parameteres of saccades detected by the Eyelink (not always reliable)
SAC = bitget(Z,ev.sacc);  

SACCADE = [];
if ~isempty(SAC)

    % beginning (wrong onset at start)
    SACCADE(:,1) = find(diff([0;SAC]) == 1);
                            
    % end (wrong offset at end)
	SACCADE(:,2) = find(diff([SAC;0]) == -1);   
                            
end

%% Smoothing for detecting saccades
% cutoff freq/NyquistFrequ
%[b,a] = butter(2,40/500,'low');
%X  = filtfilt(b, a, X);
%Y  = fil   tfilt(b, a, Y);

% Two-point central difference derivative (see Bahill et al.)
dx = NaN.*ones(length(X),1);
dy = NaN.*ones(length(Y),1);
for n = 11:(length(X)-11) % by 0.02 s
    dx(n) = (X(n+10)-X(n-10))/0.02; 
    dy(n) = (Y(n+10)-Y(n-10))/0.02; 
end

%% Expand definition of saccades by +/- 40 ms
NSAMPLE_EXPAND = 15;

ZN = fExpandSacc(Z, SACCADE, NSAMPLE_EXPAND/Ctrl.sdur,ev.blk);

dx(logical(ZN)) = NaN;
dy(logical(ZN)) = NaN;

ACCX = NaN.*ones(length(dx),1);
ACCY = NaN.*ones(length(dy),1);
for n = 11:(length(dx)-11) % by 0.02 s
    ACCX(n) = (dx(n+10)-dx(n-10))/0.02; 
    ACCY(n) = (dy(n+10)-dy(n-10))/0.02; 
end

[nsac dummy] = size(SACCADE);
if ~isempty(SACCADE)
    
    % GET X AMPLITUDE
    SACCADE(:,3) = X(SACCADE(:,2))-X(SACCADE(:,1));
    
    % GET Y AMPLITUDE
    SACCADE(:,4) = Y(SACCADE(:,2))-Y(SACCADE(:,1));
    
    % NORM: 
    SACCADE(:,5) = sqrt(SACCADE(:,3).^2 + SACCADE(:,4).^2);

    % Peak Vel
    for t = 1:nsac
        SACCADE(t,6) = max(abs(dx(SACCADE(t,1):SACCADE(t,2))));

        % Peak Vel DY 
        SACCADE(t,7) = max(abs(dy(SACCADE(t,1):SACCADE(t,2)))); 
    end
    
    % get duration
    SACCADE(:,8) = SACCADE(:,2)-SACCADE(:,1);
end

% Duration of expanded saccades:
SSAC = bitget(ZN,1);  

SSACCADE = [];
if ~isempty(SSAC)

    % beginning (wrong onset at start)
    SSACCADE(:,1) = find(diff([0;SSAC]) == 1);
                            
    % end (wrong offset at end)
	SSACCADE(:,2) = find(diff([SSAC;0]) == -1);   
                            
end         

%% get interpolated values 
y = @(x,s,i) s.*x+i;

% init
intdx = dx;
intdy = dy;

[nssac dummy] = size(SSACCADE);
for j = 1:nssac
   
   % create an interpolated trace
   if SSACCADE(1,1)>2
       
       t1 = SSACCADE(j,1);
       t2 = SSACCADE(j,2);

       y1x = dx(t1-1);
       y1y = dy(t1-1);
       
       if t2+1<recOff
       
           y2x = dx(t2+1);
           y2y = dy(t2+1);

       else
           % last non-NaN
           % if doesn't exist, last
           % value
           if isnan(dx(recOff))
               
               y2x = dx(t1-1);
               y2y = dy(t1-1);
               
           else
               
               y2x = dx(recOff);
               y2y = dy(recOff);

           end
       end

       slopex = (y2x-y1x)/(t2-t1);
       slopey = (y2y-y1y)/(t2-t1);

       % calculate intercept
       interceptx = y2x - slopex * t2;
       intercepty = y2y - slopey * t2;

       intdx((t1-1):(t2+1)) = y((t1-1):(t2+1), slopex, interceptx); 
       intdy((t1-1):(t2+1)) = y((t1-1):(t2+1), slopey, intercepty); 
       
   end
end
    