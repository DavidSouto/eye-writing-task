function ZN = sexpandsacc(Z, SACCADE, expand)
% add ones before and after a saccade
ZN = bitand(Z,1);    % saccades events

if ~isempty(SACCADE),
    for i=1:length(SACCADE(:,1)),
        j = max( [1; SACCADE(i,1)-expand]); % j cannot be zero (see previous versions)
        k = min( [length(ZN(:,1));SACCADE(i,2)+expand]);
        ZN( j:k ) = 1;
    end
end
ZN = bitor(ZN,Z);