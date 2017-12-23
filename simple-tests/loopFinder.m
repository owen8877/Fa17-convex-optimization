function signM = loopFinder(x, orig)
    [m, n] = size(x);
    path = zeros(2, m+n);
    [rowI, colI, ~] = find(x);
    rowI = [rowI; orig(1)];
    colI = [colI; orig(2)];
    
    function success = searcher(depth, u, v, searchOnRow)
        if u == orig(1) && v == orig(2) && depth > 1
            % yes, we have found the loop!
           success = true;
           return
        end
        if searchOnRow
            cList = colI(rowI == u);
            for q = cList'
                if q == v
                    continue
                end
                if searcher(depth+1, u, q, ~searchOnRow)
                    % happy!
                    path(:, depth+1) = [u; v];
                    success = true;
                    return
                end
            end
        else
            rList = rowI(colI == v);
            for q = rList'
                if q == u
                    continue
                end
                if searcher(depth+1, q, v, ~searchOnRow)
                    % happy!
                    path(:, depth+1) = [u; v];
                    success = true;
                    return
                end
            end
        end
        success = false;
    end

    searcher(0, orig(1), orig(2), false);
    path = path(:, path(1, :) > 0);
    signM = sparse(path(1, :), path(2, :), repmat([1 -1], 1, size(path, 2)/2), m, n);
end

