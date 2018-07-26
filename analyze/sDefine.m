%% That's where all parameters are defined
%clear;
format bank;

MAR     = ['s', 'o', 's', 'o',];    % Marker types
LST     = ['-';':';':';':'];        % Line Styles
MFC     = [[0 0 0];[0.5 0.5 0.5]];        % Colors

% EL2: actually frequency is detected later in sLoad, for confirmation
Ctrl.sfreq  = 1000;             % Sample Frequency in Hertz
Ctrl.sdur   = 1000/Ctrl.sfreq;  % Duration of one sample

%% screen parametesr
Ctrl.xpixel = 1280;             % number of horizontal pixels
Ctrl.ypixel = 1024;             % number of pixels vertically
Ctrl.xcm    = 39.0;             % horizontal extend of display area (cm)
Ctrl.ycm    = 29.2;             % vertical extend of display area (cm)
Ctrl.subcm  = 67.0;             % distance of subjet in cm

Ctrl.ffreq  = 85;               % frame frequency in Hertz

Ctrl.fdur   = 1000/Ctrl.ffreq;
Ctrl.xp2deg = (180/pi)*atan((Ctrl.xcm/Ctrl.xpixel) / Ctrl.subcm);
Ctrl.yp2deg = (180/pi)*atan((Ctrl.ycm/Ctrl.ypixel) / Ctrl.subcm);
Ctrl.deg2xp = 1/Ctrl.xp2deg; 

Ctrl.amplitude = 7.5;

%Ctrl.gain_crit = 0.6;  % low gain means higher or lower vertical velocity than (1-crit) a criterion too

%Ctrl.cycles_per_second = [0.2 0.15 0.1];
%Ctrl.f_frames_per_cycle = fix((1000/Ctrl.fdur)./Ctrl.cycles_per_second);
%Ctrl.critical_interval = 1000:3000; % locked to the beginning of the cycle

%Ctrl.speed = [1.0 1.2 1.4]; % internal, what does it correspond to in terms of normal walking speed?
%Ctrl.vels = [9.61 11.53 13.45]; % deg/sec

%Ctrl.range = 200;

%srcpath = 'C:\Users\Asus\Documents\src\';
addpath(genpath('../../../MATLAB/'));

%load the colorblind map
%load([srcpath 'MATLAB/FIGURES/colorblind_colormap.mat']);
Ctrl.COL = {'red','darkgreen','orange'};
% 
% for n = 1:length(COL)
%     COLMAP(n,:) = colorblind(strcmp(colornames,COL(n)),:);
% end

%% Exp 
Trl.path    = [];
Trl.exp     = 1;
Trl.trl     = 1;
Trl.num     = 1;

Trl.STR     = [];
Trl.SIDE    = [-1 1];
Trl.ANT = {'TOWARDS', 'AWAY'};

%% plotting defaults
Man.plot    = 1;
Man.sign    = [-1 1];
Man.COLOR   = ['y', 'm', 'c', 'r', 'g', 'b', 'k'];
Man.MARKER  = ['s', 'o', '*', 'x', '.', 'd', '+'];
Man.LINESTYLE = {'-';'--';':';'-.';'-';'-'; '-'};
Man.off     = 0;
Man.range   = 10; 
Man.PosRange = 10;
Man.VelRange = 20;
Man.scale_v = 0.25;
Man.datpath = '../data/';
Man.log  = [];
Man.rejected = [];

%% EYELINK 1000 EVENT TYPES
% binary
ev.msg = 4;
ev.blk = 2;
ev.sacc = 1;

%%
Err.sacc_int = 3; % saccade within critical interval
Err.blink_int = 5; % blink within critical interval
Err.msg = 6;
Err.low_gain = 7; % gain too low;
Err.num_events = 8; % incorrect number of events
Err.missing_data = 9;

Err.TYPES   = {'saccade_int' 'blink_int' 'msg' 'low gain' 'event error' 'missing data'};
Err.CODES   = [Err.sacc_int Err.blink_int Err.msg Err.low_gain Err.num_events Err.missing_data];
                 
%      fprintf(fid, '%s %d %d %d %d %d %d %2.2f\n', trl.exp, trl.sub, trl.sess, trl.num, trl.trl, trl.blk, ...
%        trl.antidx, trl.sideidx);
%    fprintf(fid, '%d %d %d %d %d %d %s %s %d\n', ... 
%                trl.exp, trl.expon, trl.sub, trl.blk, trl.num, trl.trl, trl.word, trl.word_to_write, trl.repeat);     
%    fclose(fid);

%% Conditions log-file columns 
Col.exp = 1;
Col.expon = 2;
Col.sub = 3;
Col.blk = 4;
Col.num = 5;
Col.trl = 6; 
Col.word = 7;   
Col.word_to_write = 8;
% record dwelling time ...
Col.dwelling_time = 10;
Col.control = 11;

%Col.repeat = 8;

% Letter traj (for letter and cursor traj data)
%    util.clob_name = sprintf('%sc.%ds%d.mat',util.dir_name, trl.sess, trl.trl);
%    DATA = [POS.x, POS.y, POS.t, POS.evt];
p.x =1;
p.y =2;
p.t =3;
p.evt =4;