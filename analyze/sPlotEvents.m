function sPlotEvents(Ctrl, SCR, offset, scale)
% SCR is some event
% Ctrl is a fixed parameter (or controled variable) (refresh rate, etc.)
hold on;
for i=1:length(SCR)
    plot(Ctrl.sdur*[SCR(i) SCR(i)], [(offset-scale) (offset+scale)], '-m');
end