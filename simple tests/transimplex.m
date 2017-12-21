function [X, out, val] = transimplex(~, cost, mu, nu, opts)
    [m, n] = size(cost);
    
    % generate a first feasible solution
    x = sparse(m, n);
    tmpmu = mu; tmpnu = nu;
    i = 1; j = 1;
    for l = 1:m+n-1
        if tmpmu(i) > tmpnu(j)
            x(i, j) = tmpnu(j);
            tmpmu(i) = tmpmu(i) - tmpnu(j);
            j = j + 1;
        else
            x(i, j) = tmpmu(i);
            tmpnu(j) = tmpnu(j) - tmpmu(i);
            i = i + 1;
        end
    end
    
%     function setuv(index, lastindex, dir)
%         if dir % know u, operate on v
%             indexes = find(x(index, :));
%             for op = indexes
%                 if op == lastindex
%                     continue
%                 end
%                 v(op) = cost(index, op) - u(index);
%                 setuv(op, index, false);
%             end
%         else
%             indexes = find(x(:, index))';
%             for op = indexes
%                 if op == lastindex
%                     continue
%                 end
%                 u(op) = cost(op, index) - v(index);
%                 setuv(op, index, true);
%             end
%         end
%     end
    
    itr = 1;
%     uncertainx = cell(m, 1);
%     uncertainy = cell(n, 1);
%     uncertain = java.util.HashSet(m+n);
    for i = 1:m
        uncertainx{i} = java.util.HashSet(5);
    end
    for j = 1:n
        uncertainy{j} = java.util.HashSet(5);
    end
    while true
        % compute multipliers u and v
        u = zeros(m, 1); v = zeros(1, n);
        v(end) = 0;
        
        [xii, xjj, ~] = find(x);
        undecided = sparse(xii, xjj, 1);
        
        activeList = zeros(3, m+n);
        aIndex = 0;
        [plistx, plisty, ~] = find(undecided);
        prepareList = plistx(plisty == n);
        nextSize = size(prepareList, 1);
        activeList(1, aIndex+1:aIndex+nextSize) = prepareList';
        activeList(2, aIndex+1:aIndex+nextSize) = n;
        activeList(3, aIndex+1:aIndex+nextSize) = 0;
        aIndex = aIndex + nextSize;
        
        while nnz(undecided) > 0
            tmpList = zeros(3, m+n);
            tmpIndex = 0;
            for orig = activeList(:, 1:aIndex)
                ui = orig(1);
                vi = orig(2);
                if orig(3) == 1
                    % know u and compute v
                    v(vi) = cost(ui, vi) - u(ui);
                    undecided(ui, vi) = 0;
                    [plistx, plisty, ~] = find(undecided);
                    prepareList = plistx(plisty == vi);
                    nextSize = size(prepareList, 1);
                    tmpList(1, tmpIndex+1:tmpIndex+nextSize) = prepareList';
                    tmpList(2, tmpIndex+1:tmpIndex+nextSize) = vi;
                    tmpList(3, tmpIndex+1:tmpIndex+nextSize) = 0;
                    tmpIndex = tmpIndex + nextSize;
                else
                    % know v and compute u
                    u(ui) = cost(ui, vi) - v(vi);
                    undecided(ui, vi) = 0;
                    [plistx, plisty, ~] = find(undecided);
                    prepareList = plisty(plistx == ui);
                    nextSize = size(prepareList, 1);
                    tmpList(1, tmpIndex+1:tmpIndex+nextSize) = ui;
                    tmpList(2, tmpIndex+1:tmpIndex+nextSize) = prepareList';
                    tmpList(3, tmpIndex+1:tmpIndex+nextSize) = 1;
                    tmpIndex = tmpIndex + nextSize;
                end
            end
            activeList = tmpList;
            aIndex = tmpIndex;
        end
        
        % compute relative cost
        R = cost - u - v;
        % find the most negative element
        [rowMax, colIndex] = max(-R, [], 2);
        [nrmax, nrIndex] = max(rowMax, [], 1);
        ncIndex = colIndex(nrIndex);
        if nrmax <= 0
            break
        end
        
        % assign sign to basi
        [xi, xj, ~] = find(x);
        % signM = spalloc(m, n, m+n);
        signM = sparse([xi; nrIndex], [xj; ncIndex], [zeros(size(xi, 1), 1); 1]);
        uncertain = sparse(xi, xj, 1);
        [uncertainRow, uncertainColumn] = find(uncertain);
        uindex = size(uncertainRow, 1);
%         for l = 1:size(xi, 1)
% %             if ~uncertainx.containsKey(xi(l))
% %                 uncertainx.put(xi(l), java.util.HashSet());
% %             end
%             uncertainx{xi(l)}.add(xj(l));
% %             if ~uncertainy.containsKey(xj(l))
% %                 uncertainy.put(xj(l), java.util.HashSet());
% %             end
%             uncertainy{xj(l)}.add(xi(l));
%             uncertain.add([xi(l), xj(l)]);
%         end
        
%         iterator = uncertain.iterator();
%         while ~uncertain.isEmpty()
        while nnz(uncertain) > 0
%             uItem = iterator.next();
%             ur = uItem(1);
%             uc = uItem(2);
            ur = uncertainRow(uindex);
            uc = uncertainColumn(uindex);
            
            % scan row first
%             if uncertainx{ur}.size() == 1
            if size(find(uncertain(ur, :)), 2) == 1
                % yeah, just one element left
                signM(ur, uc) = -sum(signM(ur, :));
%                 uncertainx{ur}.clear();
%                 uncertainy{uc}.remove(ur);
%                 iterator.remove();
                uncertain(ur, uc) = 0;
                [uncertainRow, uncertainColumn] = find(uncertain);
%             elseif uncertainy{uc}.size() == 1
            elseif size(find(uncertain(:, uc)), 1) == 1
                signM(ur, uc) = -sum(signM(:, uc));
%                 uncertainy{uc}.clear();
%                 uncertainx{ur}.remove(uc);
%                 iterator.remove();
                uncertain(ur, uc) = 0;
                [uncertainRow, uncertainColumn] = find(uncertain);
            end
            
%             if ~iterator.hasNext()
%                 iterator = uncertain.iterator();
%             end
            uindex = uindex - 1;
            if uindex <= 0
                uindex = size(uncertainRow, 1);
            end
        end
%         for i = 1:m
%             uncertainx{i}.clear();
%         end
%         for j = 1:n
%             uncertainy{j}.clear();
%         end
        
        % now find the min reduction
        negativePart = min(signM, 0);
        feasibleReduction = - negativePart .* x;
        [~, ~, feasibleVals] = find(feasibleReduction);
        [minReduction, ~] = min(feasibleVals);
        x = x + minReduction * signM;
        
        if mod(itr, 100) == 0
            fprintf('%d\t%.4e\n', itr, full(sum(sum(x.*cost))));
        end
        itr = itr + 1;
    end
    X = reshape(x, m*n, 1);
    out = full(sum(sum(x.*cost)));
    val = [];
end