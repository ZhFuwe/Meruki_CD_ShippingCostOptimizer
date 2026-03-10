%% 主函数
function main_Meruki_CD_ShippingCostOptimizer()
    clear; clc; close all;

    % --- 贪心算法参数 ---
    greedy_pool_size = 500;              % 生成的贪心候选解的数量
    
    % --- PSO-SA ---
    run_pso_sa_optimization = false;     % true = 深度优化, false = 快速求解
    
    % PSO-粒子群参数
    swarm_size = 50;                     % 精英粒子数量
    max_generations = 30;                % 最大演化代数
    c1 = 1.5; c2 = 1.5;                  % 学习因子
    
    % SA-模拟退火参数
    sa_iterations = 20;                  % 每个粒子进行局部搜索的迭代次数
    sa_temperature = 5000;               % 局部搜索的初始温度
    sa_alpha = 0.9;                      % 局部搜索的降温率
    
    % --- 其他参数 ---
    JPY_TO_CNY_RATE = 0.0458;            % 汇率
    max_package_count = 5;               % 最大包裹数
    alternative_threshold_ratio = 0.05;  % 备选方案阈值，相对最优成本的比例
    max_output_count = 20;               % 输出方案数量上限
    points_balance = 100;                % 积分余额
    
    % --- 检查参数 ---
    if max_package_count <= 0 || floor(max_package_count) ~= max_package_count
        error('max_package_count 必须是正整数');
    end
    if alternative_threshold_ratio < 0
        error('alternative_threshold_ratio 不能为负数');
    end
    if max_output_count <= 0 || floor(max_output_count) ~= max_output_count
        error('max_output_count 必须是正整数');
    end
    if JPY_TO_CNY_RATE <= 0
        error('JPY_TO_CNY_RATE 必须大于 0');
    end
    CNY_TO_JPY_RATE = 1 / JPY_TO_CNY_RATE; % 汇率换算


    %% --- 数据加载与预处理 ---
    try
        fprintf('正在解析数据文件...\n');
        opts = detectImportOptions('ShoppingList.csv', 'Encoding', 'UTF-8');
        opts.VariableNamingRule = 'preserve';
        numeric_cols = intersect({'序号', '价格CN', '重量', '体积', '盘数', '专辑数', '非CD金额CN', '非CD税费CN'}, opts.VariableNames);
        string_cols = intersect({'是否为单盘', '是否支持竹蜻蜓echo'}, opts.VariableNames);
        if ~isempty(numeric_cols), opts = setvartype(opts, numeric_cols, 'double'); end
        if ~isempty(string_cols), opts = setvartype(opts, string_cols, 'string'); end
        orders = readtable('ShoppingList.csv', opts);
        if ~ismember('是否支持竹蜻蜓echo', orders.Properties.VariableNames)
            orders.('是否支持竹蜻蜓echo') = repmat("Y", height(orders), 1);
        end
        if ~ismember('非CD金额CN', orders.Properties.VariableNames)
            orders.('非CD金额CN') = zeros(height(orders), 1);
        end
        if ~ismember('非CD税费CN', orders.Properties.VariableNames)
            orders.('非CD税费CN') = zeros(height(orders), 1);
        end
        orders.('非CD金额CN') = str2double(string(orders.('非CD金额CN')));
        orders.('非CD税费CN') = str2double(string(orders.('非CD税费CN')));
        orders.('非CD金额CN')(isnan(orders.('非CD金额CN'))) = 0;
        orders.('非CD税费CN')(isnan(orders.('非CD税费CN'))) = 0;
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

    %% --- 运行参数 ---
    k = max_package_count;

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
        if mod(i, 10) == 0 || i == greedy_pool_size
            current_best_cost = min(greedy_costs(1:i));
            fprintf('已生成 %d / %d 个方案，当前最优成本: %.2f JPY\n', i, greedy_pool_size, current_best_cost);
        end
    end
    [best_greedy_cost, best_greedy_idx] = min(greedy_costs);
    best_greedy_solution = greedy_solutions{best_greedy_idx};
    fprintf('贪心阶段完成。发现的最佳方案成本为: %.2f JPY\n\n', best_greedy_cost);
    
    candidate_solutions = greedy_solutions;
    candidate_costs = greedy_costs;
    
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

        candidate_solutions = [candidate_solutions; gbest_solutions(:)];
        candidate_costs = [candidate_costs; ones(length(gbest_solutions), 1) * gbest_cost];
    end

    unique_solution_bank = {};
    unique_solutions = {};
    unique_costs = [];
    for i = 1:length(candidate_solutions)
        if isinf(candidate_costs(i)) || isempty(candidate_solutions{i})
            continue;
        end
        canonical_str = getCanonicalString(candidate_solutions{i}, k, num_orders);
        existing_idx = find(strcmp(unique_solution_bank, canonical_str), 1);
        if isempty(existing_idx)
            unique_solution_bank{end+1} = canonical_str;
            unique_solutions{end+1} = candidate_solutions{i};
            unique_costs(end+1) = candidate_costs(i);
        elseif candidate_costs(i) < unique_costs(existing_idx)
            unique_solutions{existing_idx} = candidate_solutions{i};
            unique_costs(existing_idx) = candidate_costs(i);
        end
    end

    if isempty(unique_costs)
        final_solutions = {best_greedy_solution};
        final_costs = best_greedy_cost;
    else
        best_cost = min(unique_costs);
        threshold_cost = best_cost * (1 + alternative_threshold_ratio);
        keep_mask = unique_costs <= (threshold_cost + 1e-6);
        selected_solutions = unique_solutions(keep_mask);
        selected_costs = unique_costs(keep_mask);

        [sorted_selected_costs, sort_idx] = sort(selected_costs, 'ascend');
        sorted_selected_solutions = selected_solutions(sort_idx);

        output_count = min(length(sorted_selected_solutions), max_output_count);
        final_solutions = sorted_selected_solutions(1:output_count);
        final_costs = sorted_selected_costs(1:output_count);
    end

    final_cost = min(final_costs);
    threshold_cost = final_cost * (1 + alternative_threshold_ratio);
    
    %% --- 结果输出与记录 ---
    fprintf('\n--------------------------------------------------\n');
    fprintf('          *** 优化计算完成 ***\n');
    fprintf('--------------------------------------------------\n');
    [~, ~, total_points_used] = calculateCostAndDetails(final_solutions{1}, orders, k, CNY_TO_JPY_RATE, shipping_methods, weight_limits, numeric_price_data, points_balance);
    optimal_count = sum(abs(final_costs - final_cost) < 1e-6);
    alternative_count = length(final_solutions) - optimal_count;
    threshold_pct = alternative_threshold_ratio * 100;
    output_str = sprintf('查询时间: %s\n', datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    output_str = [output_str, sprintf('最终输出 %d 个方案（最优方案 %d 个，%.2f%%以内备选方案 %d 个，最多保留 %d 个）\n', length(final_solutions), optimal_count, threshold_pct, alternative_count, max_output_count)];
    output_str = [output_str, sprintf('最优成本: %.2f JPY；备选方案成本上限: %.2f JPY\n', final_cost, threshold_cost)];
    output_str = [output_str, sprintf('积分使用情况: 初始积分 %d, 使用 %d, 剩余 %d\n\n', points_balance, total_points_used, points_balance - total_points_used)];
    for i = 1:length(final_solutions)
        solution = final_solutions{i};
        solution_cost = final_costs(i);
        if final_cost > 0
            delta_pct = (solution_cost - final_cost) / final_cost * 100;
        else
            delta_pct = 0;
        end
        if abs(solution_cost - final_cost) < 1e-6
            solution_tag = '最优方案';
        else
            solution_tag = sprintf('备选方案（高 %.2f%%）', delta_pct);
        end
        [~, package_details] = calculateCostAndDetails(solution, orders, k, CNY_TO_JPY_RATE, shipping_methods, weight_limits, numeric_price_data, points_balance);
        output_str = [output_str, sprintf('=============== 方案 %d：%s ===============\n', i, solution_tag)];
        output_str = [output_str, sprintf('方案总成本: %.2f JPY\n', solution_cost)];
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
    base_packages = struct('Items', [], 'OrderCount', 0, 'TotalWeight', 0, 'TotalVolume', 0, 'TotalPriceCNY', 0, 'TotalNonCDPriceCNY', 0, 'TotalNonCDDutyCNY', 0, 'TotalDiscs', 0, 'TotalAlbums', 0, 'HasMultiDiscAlbum', false, 'SupportsEcho', true);
    base_packages = repmat(base_packages, k, 1);
    has_echo_support_col = ismember('是否支持竹蜻蜓echo', final_orders.Properties.VariableNames);
    for i = 1:length(final_solution)
        pkg_idx = final_solution(i);
        order_row = final_orders(i, :);
        base_packages(pkg_idx).Items = [base_packages(pkg_idx).Items, order_row{1, 1}];
        base_packages(pkg_idx).OrderCount = base_packages(pkg_idx).OrderCount + 1;
        base_packages(pkg_idx).TotalWeight = base_packages(pkg_idx).TotalWeight + order_row{1, 5};
        base_packages(pkg_idx).TotalVolume = base_packages(pkg_idx).TotalVolume + order_row{1, 6};
        base_packages(pkg_idx).TotalPriceCNY = base_packages(pkg_idx).TotalPriceCNY + order_row{1, 4};
        base_packages(pkg_idx).TotalNonCDPriceCNY = base_packages(pkg_idx).TotalNonCDPriceCNY + order_row{1, '非CD金额CN'};
        base_packages(pkg_idx).TotalNonCDDutyCNY = base_packages(pkg_idx).TotalNonCDDutyCNY + order_row{1, '非CD税费CN'};
        base_packages(pkg_idx).TotalDiscs = base_packages(pkg_idx).TotalDiscs + order_row{1, 7};
        base_packages(pkg_idx).TotalAlbums = base_packages(pkg_idx).TotalAlbums + order_row{1, 9};
        if strcmp(upper(strtrim(order_row{1, 8})), 'N'), base_packages(pkg_idx).HasMultiDiscAlbum = true; end
        if has_echo_support_col
            echo_support_flag = upper(strtrim(string(order_row{1, '是否支持竹蜻蜓echo'})));
            if echo_support_flag == "N", base_packages(pkg_idx).SupportsEcho = false; end
        end
    end
    opportunities = []; base_cost_sum = 0; final_package_choices = cell(k, 1);
    for i = 1:k
        if base_packages(i).OrderCount == 0 continue; end
        pkg = base_packages(i);
        [best_ineligible, best_eligible] = getBestOptions(pkg, methods, w_limits, p_data, rate);
        if isinf(best_ineligible.cost) && isinf(best_eligible.cost), total_cost=inf; package_details=[]; total_points_used=0; return; end
        if isfinite(best_ineligible.cost) && isfinite(best_eligible.cost) && best_eligible.ship > 0
            switch_cost = best_eligible.cost - best_ineligible.cost;
            if switch_cost >= 0
                opportunities(end+1,:) = [i, switch_cost, best_eligible.ship];
            end
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
    empty_package = struct('Items',[],'OrderCount',0,'TotalWeight',0,'TotalVolume',0,'TotalPriceCNY',0,'TotalNonCDPriceCNY',0,'TotalDiscs',0,'TotalAlbums',0,'HasMultiDiscAlbum',false,'BestMethod','N/A','ShippingFee',0,'CDDuty',0,'NonCDDuty',0,'Duty',0,'HandlingFee',0,'PointsUsed',0,'TotalPackageCost',0);
    package_details = repmat(empty_package, k, 1);
    for i = 1:k
        if base_packages(i).OrderCount == 0 continue; end
        final_choice = final_package_choices{i};
        package_details(i).Items=base_packages(i).Items; package_details(i).OrderCount=base_packages(i).OrderCount;
        package_details(i).TotalWeight=base_packages(i).TotalWeight; package_details(i).TotalVolume=base_packages(i).TotalVolume;
        package_details(i).TotalPriceCNY=base_packages(i).TotalPriceCNY; package_details(i).TotalNonCDPriceCNY=base_packages(i).TotalNonCDPriceCNY;
        package_details(i).TotalDiscs=base_packages(i).TotalDiscs;
        package_details(i).TotalAlbums=base_packages(i).TotalAlbums; package_details(i).HasMultiDiscAlbum=base_packages(i).HasMultiDiscAlbum;
        package_details(i).BestMethod=final_choice.method; package_details(i).ShippingFee=final_choice.ship;
        package_details(i).CDDuty=final_choice.cd_duty; package_details(i).NonCDDuty=final_choice.noncd_duty;
        package_details(i).Duty=final_choice.duty; package_details(i).HandlingFee=final_choice.handle;
        if isfield(final_choice,'points_used'), package_details(i).PointsUsed=final_choice.points_used; end
        package_details(i).TotalPackageCost=final_choice.cost;
    end
end

function [best_ineligible, best_eligible] = getBestOptions(pkg, methods, w_limits, p_data, rate)
    best_ineligible = struct('cost',inf,'method','N/A','ship',0,'duty',0,'cd_duty',0,'noncd_duty',0,'handle',0);
    best_eligible = struct('cost',inf,'method','N/A','ship',0,'duty',0,'cd_duty',0,'noncd_duty',0,'handle',0);
    for j = 1:length(methods)
        method=methods(j); 
        shipping_fee=getShippingFee(w_limits,p_data,j,pkg.TotalWeight);
        declared_price_cny = pkg.TotalPriceCNY + pkg.TotalNonCDPriceCNY;
        
        if isinf(shipping_fee) continue; end % 超重
        if strcmp(method,'竹蜻蜓Plus')&&pkg.TotalVolume>9000 continue; end % 体积限制
        if contains(method,'邮政') && declared_price_cny > 2000 continue; end % 邮政单包金额限制
        if strcmp(method,'竹蜻蜓Echo') && declared_price_cny > 2000 continue; end % 竹蜻蜓Echo单包金额限制
        if strcmp(method,'竹蜻蜓Echo') && ~pkg.SupportsEcho continue; end % 订单不支持竹蜻蜓Echo
        if strcmp(method,'竹蜻蜓Echo') && (pkg.TotalNonCDPriceCNY > 0 || pkg.TotalNonCDDutyCNY > 0) continue; end % Echo不支持非CD商品

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

        % 关税（CD税费与非CD税费分别计算）
        cd_duty_cny = 0;
        noncd_duty_cny = pkg.TotalNonCDDutyCNY;
        total_duty_cny = 0;

        if strcmp(method, '竹蜻蜓Plus')
            cd_duty_cny = 0;
            noncd_duty_cny = 0;
            total_duty_cny = 0;
        elseif strcmp(method, '竹蜻蜓') || strcmp(method, '竹蜻蜓Max')
            cd_duty_cny = pkg.TotalPriceCNY * 0.2; % 无CD免税政策
            total_duty_cny = cd_duty_cny + noncd_duty_cny;
        elseif strcmp(method, '竹蜻蜓Echo')
            if pkg.TotalDiscs > 20 % Echo适用CD免税政策
                cd_duty_cny = pkg.TotalPriceCNY * 0.2;
            end
            total_duty_cny = cd_duty_cny;
            noncd_duty_cny = 0;
        elseif contains(method, '邮政')
            if pkg.TotalDiscs > 20 % CD免税政策：20盘及以下免税
                cd_duty_cny = pkg.TotalPriceCNY * 0.2;
            end
            total_tax_cny = cd_duty_cny + noncd_duty_cny;
            if total_tax_cny > 50 % 邮政总税费门槛
                total_duty_cny = total_tax_cny;
            else
                total_duty_cny = 0;
            end
        end

        cd_duty_jpy = cd_duty_cny * rate;
        noncd_duty_jpy = noncd_duty_cny * rate;
        duty = total_duty_cny * rate;

        cost = shipping_fee+duty+handling_fee;
        current_choice = struct('cost',cost,'method',method,'ship',shipping_fee,'duty',duty,'cd_duty',cd_duty_jpy,'noncd_duty',noncd_duty_jpy,'handle',handling_fee);
        
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
        for i = 1:length(sorted_packages)
            str_parts{i} = sprintf('%d,', sorted_packages{i});
            str_parts{i}(end) = [];
        end
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
            '  - CD总价格: %.2f CNY\n' ...
            '  - 非CD总金额: %.2f CNY\n' ...
            '  - 总盘数: %d\n' ...
            '  - 专辑数: %d\n' ...
            '  - 运费(抵扣前): %.2f JPY\n' ...
            '  - 积分抵扣: -%d JPY\n' ...
            '  - 运费(抵扣后): %.2f JPY\n' ...
            '  - CD税费: %.2f JPY\n' ...
            '  - 非CD税费: %.2f JPY\n' ...
            '  - 关税: %.2f JPY\n' ...
            '  - 手续费: %d JPY\n' ...
            '  - 此包裹总计: %.2f JPY\n\n'], ...
            i, package_details(i).BestMethod, items_str, package_details(i).OrderCount, ...
            package_details(i).TotalWeight, package_details(i).TotalVolume, ...
            package_details(i).TotalPriceCNY, package_details(i).TotalNonCDPriceCNY, package_details(i).TotalDiscs, package_details(i).TotalAlbums, ...
            package_details(i).ShippingFee, ...
            package_details(i).PointsUsed, ...
            package_details(i).ShippingFee - package_details(i).PointsUsed, ...
            package_details(i).CDDuty, package_details(i).NonCDDuty, ...
            package_details(i).Duty, package_details(i).HandlingFee, package_details(i).TotalPackageCost)];
    end
end