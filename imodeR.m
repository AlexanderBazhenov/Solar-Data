function [mode, modefreq, freqs, Ss] = imodeR(S)

n = size(S,1);
bounds = repmat([1 -1], n, 1);
[Ss, idx] = sort(S(:));

bounds_sorted = bounds(idx);
freqs = cumsum(bounds_sorted);

if nargout > 2
    T = [ Ss bounds(idx) freqs ];
    [Tu, ~, Tu_ic ] = unique(T(:,[1 2]),'rows');
    count = histc(Tu_ic, 1:size(Tu,1));
    nonunique_idx = find(count > 1)';

    for k = nonunique_idx
        group_idx = find(Tu_ic == k);
        z = T(group_idx(1),2);
        mx = z*max(z*freqs(group_idx));
        T(group_idx,3) = mx;
    end

    T(:,2) = -T(:,2);
    Tu = unique(T,'rows');

    freqs = Tu(1:end-1,3);
    Ss = Tu(:,1);
end

modefreq = max(freqs);
argmax = find(freqs == modefreq);
mode = [Ss(argmax) Ss(argmax+1)];

end

