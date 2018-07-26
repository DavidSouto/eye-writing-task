% TBD: need to check event parsing works fine (sMenu)
clear;
sDefine;                                                                   % load exp info.
% Load subject data
Trl.sub     = input('sub n°: '); % subject number
Man.log = sprintf('%se%ds%d/e%ds%d.log', Man.datpath, Trl.exp, Trl.sub, Trl.exp, Trl.sub);
%DATA = textread(Man.log); %#ok<REMFF1>

fid=fopen(Man.log);
% 1 0 1 1 1 2 proved proved 0
DATA = textscan(fid,'%s','delimiter','\t');
fclose(fid);

nTrials = length(DATA{1}); % Init info
for n = 1:nTrials
    %DATA{n}=textscan(DATA{1}{n},'%d %d %d %d %d %d %s %s %d','');
    D{n}=textscan(DATA{1}{n},'%d %d %d %d %d %d %6s %6s %d','delimiter','');
    for c = 1:9 % 2 missing
        F(c).i{n} = D{n}{c};
    end
end
    
error = zeros(nTrials,1); % for recording errors

%% Load trajectory and letter info
for n = 1:nTrials                                               % trial loop
      Trl.num = n;
      sLoadLetters;
end
SDATA = dataset({[cell2mat(F(Col.sub).i'), cell2mat(F(Col.blk).i') cell2mat(F(Col.control).i'), Wo(:,1), Wo(:,2)],'sub','sess','control','error','duration'});
export(SDATA,'file',['Writing.Participant_' int2str(Trl.sub) '.txt'],'Delimiter',';');
%SDATA = dataset({[DATA(:,[Col.sub Col.sess Col.antidx]), SPAR(:,1) SPAR(:,2) dir_error error],'sub','sess','anti','lat','ampl','wdir','err'});
