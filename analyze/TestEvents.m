sDefine;

Trl.sub = 99;
Man.log = sprintf('%se%ds%d/e%ds%d.log', Man.datpath, Trl.exp, Trl.sub, Trl.exp, Trl.sub);
DATA = textread(Man.log); %#ok<REMFF1>

for num = 1:length(DATA(:,1))
    Trl.num = num;
    Man.currfile = sprintf('%se%ds%d/%db%dt%d.txt',...
    Man.datpath,Trl.exp,Trl.sub, DATA(Trl.num,Col.sess), DATA(Trl.num,Col.blk),DATA(Trl.num,Col.trl));

    if exist(Man.currfile,'file')
        TMP = textread(Man.currfile);
        %[TIME, d, d, d, X, Y, Z] 
        TIME = TMP(:,1);
        X = TMP(:,6);
        Y = TMP(:,7);
        E = TMP(:,5); % message events
        Z = TMP(:,9); % saccade events
        EVE = find(bitand(E,ev.msg));    % events 
        fprintf('EVE msg %d, unique %d, %s\n', numel(EVE), numel(unique(Z)), Man.currfile);
    end
        
    unique(Z)
end