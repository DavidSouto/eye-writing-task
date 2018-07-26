% check limits
if Trl.num > nTrials
    Trl.num = nTrials;
    msgbox('End of File');
end
trial = F(Col.trl).i{Trl.num};
session = F(Col.blk).i{Trl.num};
if trial < 1        
    trial = 1;
    msgbox('Start of File');
end

% letter file
Man.currfile = sprintf('%se%ds%d/c.%ds%d.mat',...
    Man.datpath, Trl.exp, Trl.sub, session, trial);

POS = load(Man.currfile);
VE = POS.DATA(:,p.evt);
T = POS.DATA(:,p.t);

% Find trial start and written letters
id = find(~isnan(VE) & VE~=48);
%char(VE(id))
% rel to trial start
T = diff(T(id)-T(id(1)))*1000;

% everage time per letter
written = char(F(Col.word).i{Trl.num});
target = char(F(Col.word_to_write).i{Trl.num});

if length(written)==6
    
    % word stats
    Wo(Trl.num,2) = sum(T); % duration
    Wo(Trl.num,1) = double(~strcmp(written,target)); 
    
    fprintf('Word written: %s, word to write: %s, error %d, total time %2.2f\n', written, target,Wo(Trl.num,1), Wo(Trl.num,2));

    for letter = 1:6 % to calculate error likelihood (indication of cumulated error for pre-planning)
        
        % find time first letter to next, if empty ...
        Word(Trl.num).l(letter).correct = ~strcmp(written(letter),target(letter)); 
        Word(Trl.num).l(letter).t = T(letter); % time to selection 
        c(letter) = Word(Trl.num).l(letter).correct;
        t(letter) = T(letter);
    end
    fprintf('error pattern: %d %d %d %d %d %d\n',c); 
    fprintf('time for selection: %2.0f %2.0f %2.0f %2.0f %2.0f %2.0f\n',t); 

else
    fprintf('error, word less than 6 letters\n'); %, char(F(Col.word).i{Trl.num}),char(F(Col.word_to_write).i{Trl.num}));
    W(Trl.num,1) = NaN; %~strcmp(written,target); 
    W(Trl.num,2) = NaN; %sum(T); % duration
    
    % nans
    Word(Trl.num).correct = strcmp(F(Col.word).i{Trl.num},F(Col.word_to_write).i{Trl.num}); 
    
    for letter = 1:6
        Word(Trl.num).l(letter).correct = NaN;
    end
end