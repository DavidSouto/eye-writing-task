%% Temporal envelope ..., showing the stimulus for only some 200 ms 
function [C] = fGaussEnvelope(t,t0, dur) % contrast increased by this function
    sc = dur;
    
    % define a gaussian centered on t0
    C = exp(-((t-t0)^ 2)/sc ^ 2);
    %fprintf('C %2.2f, t %d t0 %d\n',C, t, t0);
    
end
