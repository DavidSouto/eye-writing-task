function sPlotPar(Trl, Col, Err, error)
% plot error type
% if error,
%     keyboard
% end
for i=1:length(Err.CODES)
    if bitget(error,i)
        Trl.STR{length(Trl.STR)+1} = ['Error = ', Err.TYPES{i}];
    end
end

if ~isempty(Trl.STR)
    text(Trl.xstr, Trl.ystr, Trl.STR);
end