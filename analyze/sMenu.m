clear;
sDefine;

Trl.sub     = input('sub n°: '); % subject number
%subcnt = Trl.sub;
%load T; % eye movement structure

Man.plot = figure(1);
set(Man.plot, 'Position', [150, 50, 800, 650]);

Man.log = sprintf('%se%ds%d/e%ds%d.log', Man.datpath, Trl.exp, Trl.sub, Trl.exp, Trl.sub);
DATA = textread(Man.log); %#ok<REMFF1>

% Eye movement 
%MDATA = DATA(DATA(:,Col.train)==0 & DATA(:,Col.task)==0,:);        
%DATA = DATA(DATA(:,Col.train)==0 & DATA(:,Col.task)==1,:);        
error = zeros(1,length(DATA(:,1)));

%Ctrl.traj_sampling=1000;
%sGetTraj;

PlotNextTrl = uicontrol(Man.plot,'string','next Trl','Position',[205,0,90,30],'callback','Trl.num = Trl.num+1; sLoad; sPlot;');
PlotPrevTrl = uicontrol(Man.plot,'string','prev. Trl','Position',[295,0,90,30],'callback','Trl.num = Trl.num-1; sLoad; sPlot;');

