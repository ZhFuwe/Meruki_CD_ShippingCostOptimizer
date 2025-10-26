%% 主函数
function main_Meruki_CD_ShippingCostOptimizer()
    clear; clc; close all;

    % --- 贪心算法参数 ---
    greedy_pool_size = 500;          % 生成的贪心候选解的数量
    
    % --- PSO-SA ---
    run_pso_sa_optimization = false; % true = 深度优化, false = 快速求解
    
    % PSO-粒子群参数
    swarm_size = 50;                 % 精英粒子数量
    max_generations = 30;            % 最大演化代数
    c1 = 1.5; c2 = 1.5;              % 学习因子
    
    % SA-模拟退火参数
    sa_iterations = 20;              % 每个粒子进行局部搜索的迭代次数
    sa_temperature = 5000;           % 局部搜索的初始温度
    sa_alpha = 0.9;                  % 局部搜索的降温率
    
    % --- 其他参数 ---
    CNY_TO_JPY_RATE = 20.59;         % 汇率
    points_balance = 2313 ;          % 积分余额

    %% --- 数据加载与预处理 ---
    try
        fprintf('正在解析数据文件...\n');
        opts = detectImportOptions('ShoppingList.csv', 'Encoding', 'UTF-8');
        opts.VariableNamingRule = 'preserve';
        opts.VariableTypes{1} = 'double'; opts.VariableTypes{4} = 'double';
        opts.VariableTypes{5} = 'double'; opts.VariableTypes{6} = 'double';
        opts.VariableTypes{7} = 'double'; opts.VariableTypes{8} = 'string';
        opts.VariableTypes{9} = 'double';
        orders = readtable('ShoppingList.csv', opts);
        num_orders = height(orders);
        raw_price_data = readcell('ExpressPrice.csv', 'Encoding', 'UTF-8');
        shipping_methods = string(raw_price_data(1, 2:end));
        weight_limits_cell = raw_price_data(2:end, 1);
        price_grid_cell = raw_price_data(2:end, 2:end);
        valid_rows = []; weight_limits = [];
        for i = 1:length(weight_limits_cell)
            val = weight_limits_cell{i}; num_val = NaN;
            if isnumeric(val), num_val = val;
            elseif (ischar(val) || isstring(val)) && ~isempty(strtrim(val)), num_val = str2double(val); end
            if ~isnan(num_val), valid_rows = [valid_rows, i]; weight_limits = [weight_limits, num_val]; end
        end
        weight_limits = weight_limits';
        price_grid_cell = price_grid_cell(valid_rows, :);
        numeric_price_data = ones(size(price_grid_cell)) * inf;
        for i = 1:size(price_grid_cell, 1)
            for j = 1:size(price_grid_cell, 2)
                val = price_grid_cell{i, j};
                if isnumeric(val), numeric_price_data(i, j) = val;
                elseif (ischar(val) || isstring(val)) && ~strcmpi(val, 'E')
                    num_val = str2double(val);
                    if ~isnan(num_val), numeric_price_data(i, j) = num_val; end
                end
            end
        end
        fprintf('数据解析完成。\n\n');
    catch ME
        disp('错误：数据文件加载或解析失败。'); disp(ME.message); return;
    end

    %% --- 用户输入 ---
    prompt = {'请输入您希望合并包裹的最大数量:'};
    dlgtitle = '设置'; dims = [1 50]; definput = {'3'};
    answer = inputdlg(prompt, dlgtitle, dims, definput);
    if isempty(answer), disp('用户取消了操作。'); return; end
    k = str2double(answer{1});
    if isnan(k) || k <= 0 || floor(k) ~= k, disp('错误：请输入一个正整数。'); return; end

    %% --- 算法主流程 ---
    fprintf('--- 阶段一：贪心策略生成 %d 个方案 ---\n', greedy_pool_size);
    greedy_strategies = {'descend weight', 'descend volume', 'descend price', 'ascend weight', 'ascend volume'};
    greedy_solutions = cell(greedy_pool_size, 1);
    greedy_costs = ones(greedy_pool_size, 1) * inf;
    for i = 1:greedy_pool_size
        strategy = greedy_strategies{randi(length(greedy_strategies))};
        [solution, cost] = createGreedyInitialSolution(orders, k, CNY_TO_JPY_RATE, shipping_methods, weight_limits, numeric_price_data, strategy, points_balance);
        greedy_solutions{i} = solution;
        greedy_costs(i) = cost;
        if mod(i, 100) == 0, fprintf('已生成 %d / %d 个方案...\n', i, greedy_pool_size); end
    end
    [best_greedy_cost, best_greedy_idx] = min(greedy_costs);
    best_greedy_solution = greedy_solutions{best_greedy_idx};
    fprintf('贪心阶段完成。发现的最佳方案成本为: %.2f JPY\n\n', best_greedy_cost);
    
    final_solutions = {best_greedy_solution};
    final_cost = best_greedy_cost;
    
    if run_pso_sa_optimization
        fprintf('--- 阶段二：PSO-SA 深度优化 ---\n');
        [sorted_costs, sorted_indices] = sort(greedy_costs);
        sorted_solutions = greedy_solutions(sorted_indices);
        elite_particles = struct('position', [], 'cost', inf, 'pbest_position', [], 'pbest_cost', inf);
        elite_count = 0;
        solution_bank = {};
        for i = 1:greedy_pool_size
            if elite_count >= swarm_size, break; end
            canonical_str = getCanonicalString(sorted_solutions{i}, k, num_orders);
            if ~ismember(canonical_str, solution_bank)
                elite_count = elite_count + 1;
                elite_particles(elite_count).position = sorted_solutions{i};
                elite_particles(elite_count).cost = sorted_costs(i);
                elite_particles(elite_count).pbest_position = sorted_solutions{i};
                elite_particles(elite_count).pbest_cost = sorted_costs(i);
                solution_bank{end+1} = canonical_str;
            end
        end
        fprintf('已筛选 %d 个精英粒子。\n\n', elite_count);
        particles = elite_particles;
        gbest_cost = particles(1).pbest_cost;
        gbest_solutions = {particles(1).pbest_position};
        gbest_bank = {getCanonicalString(particles(1).pbest_position, k, num_orders)};
        gbest_history = [gbest_cost];
        for gen = 1:max_generations
            for i = 1:elite_count
                particle = particles(i);
                new_position = particle.position;
                for order_idx = 1:num_orders
                    if rand() < c1*rand() && new_position(order_idx) ~= particle.pbest_position(order_idx), new_position(order_idx) = particle.pbest_position(order_idx); end
                    if rand() < c2*rand() && ~isempty(gbest_solutions) && new_position(order_idx) ~= gbest_solutions{1}(order_idx), new_position(order_idx) = gbest_solutions{1}(order_idx); end
                end
                sa_current_solution = new_position;
                T = sa_temperature;
                for sa_iter = 1:sa_iterations
                    sa_new_solution = generateAggressiveNeighbor(sa_current_solution, k, orders);
                    new_cost = calculateCostAndDetails(sa_new_solution, orders, k, CNY_TO_JPY_RATE, shipping_methods, weight_limits, numeric_price_data, points_balance);
                    current_cost = calculateCostAndDetails(sa_current_solution, orders, k, CNY_TO_JPY_RATE, shipping_methods, weight_limits, numeric_price_data, points_balance);
                    delta_cost = new_cost - current_cost;
                    if delta_cost < 0 || rand() < exp(-delta_cost / T), sa_current_solution = sa_new_solution; end
                    T = T * sa_alpha;
                end
                particles(i).position = sa_current_solution;
                particles(i).cost = calculateCostAndDetails(sa_current_solution, orders, k, CNY_TO_JPY_RATE, shipping_methods, weight_limits, numeric_price_data, points_balance);
                if particles(i).cost < particles(i).pbest_cost, particles(i).pbest_position = particles(i).position; particles(i).pbest_cost = particles(i).cost; end
                if particles(i).pbest_cost < gbest_cost
                    gbest_cost = particles(i).pbest_cost;
                    gbest_solutions = {particles(i).pbest_position};
                    gbest_bank = {getCanonicalString(particles(i).pbest_position, k, num_orders)};
                elseif abs(particles(i).pbest_cost - gbest_cost) < 1e-6
                    canonical_str = getCanonicalString(particles(i).pbest_position, k, num_orders);
                    if ~ismember(canonical_str, gbest_bank), gbest_solutions{end+1} = particles(i).pbest_position; gbest_bank{end+1} = canonical_str; end
                end
            end
            gbest_history = [gbest_history, gbest_cost];
            fprintf('第 %d 代演化完成, 当前最优总成本: %.2f JPY, 已发现 %d 个最优方案\n', gen, gbest_cost, length(gbest_solutions));
        end
    end
    
    %% --- 结果输出与记录 ---
    fprintf('\n--------------------------------------------------\n');
    fprintf('          *** 优化计算完成 ***\n');
    fprintf('--------------------------------------------------\n');
    [~, ~, total_points_used] = calculateCostAndDetails(final_solutions{1}, orders, k, CNY_TO_JPY_RATE, shipping_methods, weight_limits, numeric_price_data, points_balance);
    output_str = sprintf('查询时间: %s\n', datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    output_str = [output_str, sprintf('最终找到 %d 个成本为 %.2f JPY 的最优方案\n', length(final_solutions), final_cost)];
    output_str = [output_str, sprintf('积分使用情况: 初始积分 %d, 使用 %d, 剩余 %d\n\n', points_balance, total_points_used, points_balance - total_points_used)];
    for i = 1:length(final_solutions)
        solution = final_solutions{i};
        [~, package_details] = calculateCostAndDetails(solution, orders, k, CNY_TO_JPY_RATE, shipping_methods, weight_limits, numeric_price_data, points_balance);
        output_str = [output_str, sprintf('=============== 方案 %d ===============\n', i)];
        output_str = [output_str, formatOutput(package_details, k)];
    end
    disp(output_str);
    log_filename = sprintf('%s-Cost_%.2f.txt', datestr(now, 'yyyy-mm-dd-HH-MM'), final_cost);
    try
        fileID = fopen(log_filename, 'w', 'n', 'UTF-8'); fprintf(fileID, '%s', output_str); fclose(fileID);
        fprintf('结果已成功保存到新文件: %s\n', log_filename);
    catch ME
        fprintf('错误：无法写入日志文件: %s\n', log_filename); disp(ME.message);
    end
end

%% --- 辅助函数 ---

function [initial_solution, initial_cost] = createGreedyInitialSolution(orders, k, rate, methods, w_limits, p_data, strategy, points)
    num_orders = height(orders); initial_solution = zeros(1, num_orders);
    parts = strsplit(strategy); direction = parts{1}; property = parts{2};
    col_idx = 5; if strcmp(property, 'volume'), col_idx = 6; elseif strcmp(property, 'price'), col_idx = 4; end
    [~, sorted_indices] = sort(orders{:, col_idx}, direction);
    sorted_indices = sorted_indices(randperm(num_orders));
    for i = 1:num_orders
        order_idx = sorted_indices(i); best_package_idx = -1; min_total_cost = inf;
        for j = 1:k
            temp_solution = initial_solution; temp_solution(order_idx) = j;
            cost = calculateCostAndDetails(temp_solution, orders, k, rate, methods, w_limits, p_data, points);
            if cost < min_total_cost, min_total_cost = cost; best_package_idx = j; end
        end
        initial_solution(order_idx) = best_package_idx;
    end
    initial_cost = calculateCostAndDetails(initial_solution, orders, k, rate, methods, w_limits, p_data, points);
end

function [total_cost, package_details, total_points_used] = calculateCostAndDetails(solution, orders, k, rate, methods, w_limits, p_data, points)
    if isempty(solution) || all(solution == 0), total_cost = 0; package_details = []; total_points_used = 0; return; end
    final_solution = solution; final_orders = orders;
    if any(solution == 0)
        active_indices = find(solution ~= 0);
        final_solution = solution(active_indices);
        final_orders = orders(active_indices, :);
    end
    base_packages = struct('Items', [], 'OrderCount', 0, 'TotalWeight', 0, 'TotalVolume', 0, 'TotalPriceCNY', 0, 'TotalDiscs', 0, 'TotalAlbums', 0, 'HasMultiDiscAlbum', false);
    base_packages = repmat(base_packages, k, 1);
    for i = 1:length(final_solution)
        pkg_idx = final_solution(i);
        order_row = final_orders(i, :);
        base_packages(pkg_idx).Items = [base_packages(pkg_idx).Items, order_row{1, 1}];
        base_packages(pkg_idx).OrderCount = base_packages(pkg_idx).OrderCount + 1;
        base_packages(pkg_idx).TotalWeight = base_packages(pkg_idx).TotalWeight + order_row{1, 5};
        base_packages(pkg_idx).TotalVolume = base_packages(pkg_idx).TotalVolume + order_row{1, 6};
        base_packages(pkg_idx).TotalPriceCNY = base_packages(pkg_idx).TotalPriceCNY + order_row{1, 4};
        base_packages(pkg_idx).TotalDiscs = base_packages(pkg_idx).TotalDiscs + order_row{1, 7};
        base_packages(pkg_idx).TotalAlbums = base_packages(pkg_idx).TotalAlbums + order_row{1, 9};
        if strcmp(upper(strtrim(order_row{1, 8})), 'N'), base_packages(pkg_idx).HasMultiDiscAlbum = true; end
    end
    opportunities = []; base_cost_sum = 0; final_package_choices = cell(k, 1);
    for i = 1:k
        if base_packages(i).OrderCount == 0 continue; end
        pkg = base_packages(i);
        [best_ineligible, best_eligible] = getBestOptions(pkg, methods, w_limits, p_data, rate);
        if isinf(best_ineligible.cost) && isinf(best_eligible.cost), total_cost=inf; package_details=[]; total_points_used=0; return; end
        if ~isinf(best_eligible.cost) && best_eligible.ship > 0
            switch_cost = best_eligible.cost - best_ineligible.cost;
            opportunities(end+1,:) = [i, switch_cost, best_eligible.ship];
        end
        if best_ineligible.cost <= best_eligible.cost, final_package_choices{i} = best_ineligible; base_cost_sum = base_cost_sum + best_ineligible.cost;
        else, final_package_choices{i} = best_eligible; base_cost_sum = base_cost_sum + best_eligible.cost; end
    end
    remaining_points = points;
    if ~isempty(opportunities)
        opportunities = sortrows(opportunities, 2);
        for opp_idx = 1:size(opportunities,1)
            if remaining_points <= 0 break; end
            pkg_idx = opportunities(opp_idx, 1);
            switch_cost = opportunities(opp_idx, 2);
            points_needed = opportunities(opp_idx, 3);
            points_to_use = min(remaining_points, points_needed);
            if switch_cost < points_to_use
                base_cost_sum = base_cost_sum - (points_to_use - switch_cost);
                remaining_points = remaining_points - points_to_use;
                [~, best_eligible] = getBestOptions(base_packages(pkg_idx), methods, w_limits, p_data, rate);
                best_eligible.cost = best_eligible.cost - points_to_use;
                best_eligible.points_used = points_to_use;
                final_package_choices{pkg_idx} = best_eligible;
            end
        end
    end
    total_cost = base_cost_sum;
    total_points_used = points - remaining_points;
    empty_package = struct('Items',[],'OrderCount',0,'TotalWeight',0,'TotalVolume',0,'TotalPriceCNY',0,'TotalDiscs',0,'TotalAlbums',0,'HasMultiDiscAlbum',false,'BestMethod','N/A','ShippingFee',0,'Duty',0,'HandlingFee',0,'PointsUsed',0,'TotalPackageCost',0);
    package_details = repmat(empty_package, k, 1);
    for i = 1:k
        if base_packages(i).OrderCount == 0 continue; end
        final_choice = final_package_choices{i};
        package_details(i).Items=base_packages(i).Items; package_details(i).OrderCount=base_packages(i).OrderCount;
        package_details(i).TotalWeight=base_packages(i).TotalWeight; package_details(i).TotalVolume=base_packages(i).TotalVolume;
        package_details(i).TotalPriceCNY=base_packages(i).TotalPriceCNY; package_details(i).TotalDiscs=base_packages(i).TotalDiscs;
        package_details(i).TotalAlbums=base_packages(i).TotalAlbums; package_details(i).HasMultiDiscAlbum=base_packages(i).HasMultiDiscAlbum;
        package_details(i).BestMethod=final_choice.method; package_details(i).ShippingFee=final_choice.ship;
        package_details(i).Duty=final_choice.duty; package_details(i).HandlingFee=final_choice.handle;
        if isfield(final_choice,'points_used'), package_details(i).PointsUsed=final_choice.points_used; end
        package_details(i).TotalPackageCost=final_choice.cost;
    end
end

function [best_ineligible, best_eligible] = getBestOptions(pkg, methods, w_limits, p_data, rate)
    best_ineligible = struct('cost',inf,'method','N/A','ship',0,'duty',0,'handle',0);
    best_eligible = struct('cost',inf,'method','N/A','ship',0,'duty',0,'handle',0);
    for j = 1:length(methods)
        method=methods(j); 
        shipping_fee=getShippingFee(w_limits,p_data,j,pkg.TotalWeight);
        
        if isinf(shipping_fee) continue; end % 超重
        if strcmp(method,'竹蜻蜓Plus')&&pkg.TotalVolume>9000 continue; end % 体积限制

        % 订单数量限制
        if (strcmp(method, '竹蜻蜓Max') || strcmp(method, '竹蜻蜓Echo')) && pkg.OrderCount > 14
            continue;
        end
        
        % 手续费
        duty=0; 
        handling_fee=0;
        weight_g = pkg.TotalWeight;

        if strcmp(method, '邮政-航空件') || strcmp(method, '邮政-海运')
            if weight_g > 2000 && weight_g <= 5000 % 2-5kg
                handling_fee = 200;
            end
        elseif strcmp(method, '竹蜻蜓Max') || strcmp(method, '竹蜻蜓Echo')
            if weight_g <= 3000
                handling_fee = 150;
            elseif weight_g > 3000 && weight_g <= 5000
                handling_fee = 200;
            end
        elseif strcmp(method, '邮政-可追踪小航空')
            handling_fee = 150;
        elseif strcmp(method, '竹蜻蜓Plus')
            handling_fee = 200;
        end

        % 关税
        if contains(method,'邮政')
            if pkg.HasMultiDiscAlbum, if pkg.TotalAlbums>3, duty=pkg.TotalPriceCNY*0.2*rate; end
            else, if pkg.TotalDiscs>20, duty=pkg.TotalPriceCNY*0.2*rate; end
            end
        elseif contains(method,'竹蜻蜓')
            if ~strcmp(method,'竹蜻蜓Echo') % 关税 for 竹蜻蜓, Plus, Max
                if pkg.TotalPriceCNY*0.2>50, duty=pkg.TotalPriceCNY*0.2*rate; end
            else % 关税 for 竹蜻蜓Echo
                if pkg.TotalDiscs>20&&pkg.TotalDiscs>0
                    avg_p=pkg.TotalPriceCNY/pkg.TotalDiscs; 
                    tax_p=(pkg.TotalDiscs-20)*avg_p; 
                    pot_d=tax_p*0.2; 
                    if pot_d>=50, duty=pot_d*rate; end
                end
            end
        end

        cost = shipping_fee+duty+handling_fee;
        current_choice = struct('cost',cost,'method',method,'ship',shipping_fee,'duty',duty,'handle',handling_fee);
        
        if contains(method,'竹蜻蜓'), if cost<best_ineligible.cost, best_ineligible=current_choice; end
        else, if cost<best_eligible.cost, best_eligible=current_choice; end; end
    end
end

function canonical_str = getCanonicalString(solution, k, num_orders)
    packages = cell(1, k);
    for i = 1:num_orders, pkg_idx = solution(i); packages{pkg_idx} = [packages{pkg_idx}, i]; end
    non_empty_packages = {};
    for i = 1:k, if ~isempty(packages{i}), non_empty_packages{end+1} = sort(packages{i}); end; end
    if ~isempty(non_empty_packages)
        first_elements = cellfun(@(c) c(1), non_empty_packages);
        [~, sort_indices] = sort(first_elements);
        sorted_packages = non_empty_packages(sort_indices);
        str_parts = cell(1, length(sorted_packages));
        for i = 1:length(sorted_packages), str_parts{i} = strjoin(string(sorted_packages{i}), ','); end
        canonical_str = strjoin(str_parts, '_');
    else, canonical_str = ''; end
end

function fee = getShippingFee(weight_limits, numeric_price_data, method_idx, weight_g)
    weight_kg = weight_g / 1000;
    idx_row = find(weight_limits >= weight_kg, 1, 'first');
    if isempty(idx_row), fee = inf;
    else, fee = numeric_price_data(idx_row, method_idx); end
end

function new_solution = generateAggressiveNeighbor(solution, k, orders)
    new_solution = solution; num_orders = length(solution); op_rand = rand();
    if op_rand < 0.3 && k > 1, pkg_indices=randperm(k,2); o1s=find(solution==pkg_indices(1)); o2s=find(solution==pkg_indices(2)); if ~isempty(o1s)&&~isempty(o2s), o1=o1s(randi(length(o1s))); o2=o2s(randi(length(o2s))); new_solution(o1)=pkg_indices(2); new_solution(o2)=pkg_indices(1); end
    elseif op_rand < 0.6 && k > 1, [~,p_max]=findPackageByProperty(solution,orders,k,'weight','max'); [~,p_min]=findPackageByProperty(solution,orders,k,'weight','min'); if p_max~=p_min, o_move=find(solution==p_max); if length(o_move)>1, o_idx=o_move(randi(length(o_move))); new_solution(o_idx)=p_min; end; end
    elseif op_rand < 0.8 && k > 2, p=randperm(k,3); o1s=find(solution==p(1)); o2s=find(solution==p(2)); o3s=find(solution==p(3)); if ~isempty(o1s)&&~isempty(o2s)&&~isempty(o3s), o1=o1s(randi(length(o1s))); o2=o2s(randi(length(o2s))); o3=o3s(randi(length(o3s))); new_solution(o1)=p(2); new_solution(o2)=p(3); new_solution(o3)=p(1); end
    else, o_move=randi(num_orders); curr_p=new_solution(o_move); new_p=randi(k); while new_p==curr_p&&k>1, new_p=randi(k); end; new_solution(o_move)=new_p; end
end

function [val, idx] = findPackageByProperty(solution, orders, k, property, find_type)
    vals = zeros(1, k); col_idx = 5; if strcmp(property,'volume'),col_idx=6;end; pkgs=unique(solution);
    for i=1:k, if ismember(i,pkgs), vals(i)=sum(orders{solution==i,col_idx}); else, if strcmp(find_type,'max'),vals(i)=-inf;else,vals(i)=inf;end;end;end
    if strcmp(find_type,'max'),[val,idx]=max(vals);else,[val,idx]=min(vals);end
end

function output_str = formatOutput(package_details, k)
    output_str = '';
    for i=1:k
        if isempty(package_details(i).Items) continue; end
        items_str = strjoin(string(sort(unique(package_details(i).Items))), ' ');
        output_str=[output_str, sprintf(['--- 包裹 %d ---\n' ...
            '  - 最佳物流方式: %s\n' ...
            '  - 包裹内含订单序号: %s\n' ...
            '  - 订单数量: %d\n' ...
            '  - 总重量: %.2f g\n' ...
            '  - 总体积: %.2f\n' ...
            '  - 总价格: %.2f CNY\n' ...
            '  - 总盘数: %d\n' ...
            '  - 专辑数: %d\n' ...
            '  - 运费(抵扣前): %.2f JPY\n' ...
            '  - 积分抵扣: -%d JPY\n' ...
            '  - 运费(抵扣后): %.2f JPY\n' ...
            '  - 关税: %.2f JPY\n' ...
            '  - 手续费: %d JPY\n' ...
            '  - 此包裹总计: %.2f JPY\n\n'], ...
            i, package_details(i).BestMethod, items_str, package_details(i).OrderCount, ...
            package_details(i).TotalWeight, package_details(i).TotalVolume, ...
            package_details(i).TotalPriceCNY, package_details(i).TotalDiscs, package_details(i).TotalAlbums, ...
            package_details(i).ShippingFee, ...
            package_details(i).PointsUsed, ...
            package_details(i).ShippingFee - package_details(i).PointsUsed, ...
            package_details(i).Duty, package_details(i).HandlingFee, package_details(i).TotalPackageCost)];
    end
end