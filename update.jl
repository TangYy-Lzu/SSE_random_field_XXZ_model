module SSEUpdate
export diagonalupdate, loopupdate, adjustcutoff

# 对角更新、循环更新以及调整算符串长度的函数
# 在每个函数中增加了详细注释，解释每个步骤的意义和参数用法

function diagonalupdate(nh::Int64, mm::Int64, opstring::Vector{Int64}, spin::Vector{Int64}, spin_mm::Matrix{Int64}, bsites::Matrix{Int64}, nb::Int64, beta::Float64, Δ::Float64, J::Float64, h::Vector{Float64}, W::Float64)
    # nh: 当前对角算符数目, mm: 算符字符串长度, opstring: 存放算符的数组
    # spin: 自旋数组, bsites: 键对应的自旋编号, nb: 键的总数, beta: 逆温
    # 根据概率决定在空位添加或移除对角算符，并在非对角算符处进行自旋翻转
    upprob::Float64 = (W / 2 + abs(Δ / 4) - Δ / 4) * beta * nb  # 添加对角操作的概率因子，对角操作贡献一个1/2的因子
    downprob::Float64 = (W / 2 + abs(Δ / 4) - Δ / 4) * beta * nb
    difprob::Float64 = (W / 2 + abs(Δ / 4) + Δ / 4) * beta * nb
    """
    由于我们对所有bond进行操作，因此对角算符的添加和移除并不是以相同概率进行的，
    对对角算符的添加次数等于(mm-nh)，对对角算符的移除次数(nh+1)，也就是说添加对角算符的概率为(mm-nh)/(mm-nh+nh+1)，
    而选择一个空算符移除的概率正好是1/(mm-nh)，添加一个非对角算符的概率也是(1/nh+1)，总的概率正好一样，因此不用考虑我们选择在哪个字符串添加
    直接计算我们选好了一个位置，开始计算细致平衡就行。
    """
    @inbounds for i = 1:mm
        op::Int64 = opstring[i]
        # add diagonal operator
        if op == 0
            # 空位置：尝试添加对角算符 ，在某个δt上添加一个bond
            b = Int(trunc(nb * rand())) + 1  # 随机选择一个键， trunc(x) 是 Julia 的内置函数，用于对浮点数 x 进行截断（取整数部分），即舍弃小数部分。
            if spin[bsites[1, b]] != spin[bsites[2, b]]# 故意让整个哈密顿量加上1/4，放入对角算符中以让它的作用为0
                if spin[bsites[1, b]] == 1 # 上下自旋
                    p_acc = (difprob + (h[bsites[1, b]] / 4 - h[bsites[2, b]] / 4) * beta * nb) / (mm - nh)
                    if (p_acc >= 1.0 || p_acc >= rand())
                        opstring[i] = 2 * b  # 记录添加了对角算符（2*b 偶数表示对角操作符，b代表作用在哪个bond）
                        spin_mm[1, i] = 1  # 记录自旋状态
                        spin_mm[2, i] = -1  # 记录自旋状态
                        spin_mm[3, i] = 1  # 记录自旋状态
                        spin_mm[4, i] = -1  # 记录自旋状态
                        nh = nh + 1         # 增加当前对角算符数目
                    end
                else
                    p_acc = (difprob - (h[bsites[1, b]] / 4 - h[bsites[2, b]] / 4) * beta * nb) / (mm - nh) # 下上自旋
                    if (p_acc >= 1.0 || p_acc >= rand())
                        opstring[i] = 2 * b  # 记录添加了对角算符（2*b 偶数表示对角操作符，b代表作用在哪个bond）
                        spin_mm[1, i] = -1  # 记录自旋状态
                        spin_mm[2, i] = 1  # 记录自旋状态
                        spin_mm[3, i] = -1  # 记录自旋状态
                        spin_mm[4, i] = 1  # 记录自旋状态
                        nh = nh + 1         # 增加当前对角算符数目
                    end
                end
            elseif spin[bsites[1, b]] == 1
                # 当两个自旋同为向上
                p_acc = (upprob + (h[bsites[1, b]] / 4 + h[bsites[2, b]] / 4) * beta * nb) / (mm - nh)
                if (p_acc >= 1.0 || p_acc >= rand())
                    opstring[i] = 2 * b  # 记录添加了对角算符（2*b 偶数表示对角操作符，b代表作用在哪个bond）向上
                    spin_mm[1, i] = 1  # 记录自旋状态
                    spin_mm[2, i] = 1  # 记录自旋状态
                    spin_mm[3, i] = 1  # 记录自旋状态
                    spin_mm[4, i] = 1  # 记录自旋状态
                    nh = nh + 1         # 增加当前对角算符数目
                end
            else
                # 当两个自旋同为向下
                p_acc = (downprob - (h[bsites[1, b]] / 4 + h[bsites[2, b]] / 4) * beta * nb) / (mm - nh)
                if (p_acc >= 1.0 || p_acc >= rand())
                    opstring[i] = 2 * b  # 记录添加了对角算符（2*b 偶数表示对角操作符，b代表作用在哪个bond）
                    spin_mm[1, i] = -1  # 记录自旋状态
                    spin_mm[2, i] = -1  # 记录自旋状态
                    spin_mm[3, i] = -1  # 记录自旋状态
                    spin_mm[4, i] = -1  # 记录自旋状态
                    nh = nh + 1         # 增加当前对角算符数目
                end
            end
            # remove diagonal operator
        elseif mod(op, 2) == 0
            b = op / 2
            site_1 = bsites[1, b]
            site_2 = bsites[2, b]
            # 当前为对角算符：尝试移除已存在的对角操作
            if spin[site_1] != spin[site_2]
                if spin[bsites[1, b]] == 1 # 上下自旋
                    p_acc = (difprob + (h[bsites[1, b]] / 4 - h[bsites[2, b]] / 4) * beta * nb)
                    p_ref = (mm - nh + 1) / p_acc
                    if (p_ref >= 1.0 || p_ref >= rand())
                        opstring[i] = 0    # 移除操作 0 为单位算符
                        spin_mm[1, i] = 1  # 记录自旋状态
                        spin_mm[2, i] = -1  # 记录自旋状态
                        spin_mm[3, i] = 1  # 记录自旋状态
                        spin_mm[4, i] = -1  # 记录自旋状态
                        nh = nh - 1        # 更新对角算符数目
                    end
                else # 下上自旋
                    p_acc = (difprob - (h[bsites[1, b]] / 4 - h[bsites[2, b]] / 4) * beta * nb) # 下上自旋
                    p_ref = (mm - nh + 1) / p_acc
                    if (p_ref >= 1.0 || p_ref >= rand())
                        opstring[i] = 0    # 移除操作 0 为单位算符
                        spin_mm[1, i] = -1  # 记录自旋状态
                        spin_mm[2, i] = 1  # 记录自旋状态
                        spin_mm[3, i] = -1  # 记录自旋状态
                        spin_mm[4, i] = 1  # 记录自旋状态
                        nh = nh - 1        # 更新对角算符数目
                    end
                end
            elseif spin[bsites[1, b]] == 1
                # 当两个自旋同为向上
                p_acc = (upprob + (h[bsites[1, b]] / 4 + h[bsites[2, b]] / 4) * beta * nb)
                p_ref = (mm - nh + 1) / p_acc
                if (p_ref >= 1.0 || p_ref >= rand())
                    opstring[i] = 0    # 移除操作 0 为单位算符
                    spin_mm[1, i] = 1  # 记录自旋状态
                    spin_mm[2, i] = 1  # 记录自旋状态
                    spin_mm[3, i] = 1  # 记录自旋状态
                    spin_mm[4, i] = 1  # 记录自旋状态
                    nh = nh - 1        # 更新对角算符数目
                end
            else
                # 当两个自旋同为向下
                p_acc = (downprob - (h[bsites[1, b]] / 4 + h[bsites[2, b]] / 4) * beta * nb)
                p_ref = (mm - nh + 1) / p_acc
                if (p_ref >= 1.0 || p_ref >= rand())
                    opstring[i] = 0    # 移除操作 0 为单位算符
                    spin_mm[1, i] = -1  # 记录自旋状态
                    spin_mm[2, i] = -1  # 记录自旋状态
                    spin_mm[3, i] = -1  # 记录自旋状态
                    spin_mm[4, i] = -1  # 记录自旋状态
                    nh = nh - 1        # 更新对角算符数目
                end
            end
        else
            # 非对角算符：执行自旋翻转操作，翻转后保存这个状态，这是为了得到链表下一时刻的自旋，用来判断是否可以放入对角算符
            # 所有都执行完以后，应该还是和原来一样
            b::Int64 = fld(op, 2)
            if spin[bsites[1, b]] == 1
                # 上下自旋
                spin_mm[1, i] = 1  # 记录自旋状态
                spin_mm[2, i] = -1  # 记录自旋状态
                spin_mm[3, i] = -1  # 记录自旋状态
                spin_mm[4, i] = 1  # 记录自旋状态
                spin[bsites[1, b]] = -spin[bsites[1, b]]
                spin[bsites[2, b]] = -spin[bsites[2, b]]
            else # 下上自旋
                spin_mm[1, i] = -1  # 记录自旋状态
                spin_mm[2, i] = 1  # 记录自旋状态
                spin_mm[3, i] = 1  # 记录自旋状态
                spin_mm[4, i] = -1  # 记录自旋状态
                spin[bsites[1, b]] = -spin[bsites[1, b]]
                spin[bsites[2, b]] = -spin[bsites[2, b]]
            end
        end
    end
    return nh
end

function loopupdate(mm::Int64, nn::Int64, opstring::Vector{Int64}, bsites::Matrix{Int64}, vertexlist::Vector{Int64}, spin::Vector{Int64}, spin_mm::Matrix{Int64}, Δ::Float64, J::Float64, h::Vector{Float64}, W::Float64)
    # 初始化：firstspinop 用于记录每个自旋第一次出现的顶点索引,
    # lastspinop 用于记录每个自旋最后一次出现的顶点索引
    firstspinop::Vector{Int64} = fill(-1, nn)
    lastspinop::Vector{Int64} = fill(-1, nn)

    # 修改自 sandvik ssebasic.f90 代码
    # 主要修改：顶点腿标号从 1,2,3,4 而非 0,1,2,3

    # ============================================================
    # 构建链式顶点列表
    # ============================================================
    @inbounds for v0 in range(start=1, step=4, stop=4 * mm - 3) # 遍历所有顶点编号 v0，每次以步长为 4 递增，用于处理随机级数展开蒙特卡洛（SSE）方法中每个算符对应的顶点组。   
        """
        第 1 个算符对应顶点编号 1, 2, 3, 4。
        第 2 个算符对应顶点编号 5, 6, 7, 8。
        第 3 个算符对应顶点编号 9, 10, 11, 12。
        """
        op = opstring[cld(v0, 4)] # cld(v0, 4) 是一个自定义函数，表示向上取整到 4 的倍数
        if op != 0  # 如果该位置有算符
            # 获取当前算符对应的 bond 编号
            b = fld(op, 2) # 向下取整，偶数为对角算符，奇数非对角，模表示位置，即编号
            # 从 bsites 矩阵中获得该 bond 对应的两个自旋索引，b位置的键作用的两个自旋的位置为s1,s2
            s1 = bsites[1, b]
            s2 = bsites[2, b]
            # 获取当前自旋 s1 和 s2 最近一次出现顶点索引
            v1 = lastspinop[s1]
            v2 = lastspinop[s2]
            # 如果 s1 已有上一个顶点（已经和前面的顶点链接），则建立链式链接
            if v1 != -1 # 表示自旋s1已经出现过，并且v1代表上次出现与v0作用在相同自旋的顶点的编号
                vertexlist[v1] = v0 # 链接v0与v1
                vertexlist[v0] = v1 # 相互链接
            else
                # s1 第一次出现，把当前顶点记录下来，firstspinop用于记录每个自旋第一次出现的顶点索引
                firstspinop[s1] = v0
            end
            # 对 s2 做类似处理；注意顶点标号偏移
            if v2 != -1
                vertexlist[v2] = v0 + 1
                vertexlist[v0+1] = v2
            else
                firstspinop[s2] = v0 + 1
            end
            # 更新 s1 和 s2 的最后出现顶点索引，v0和v0+4是当前顶点，作用上算符后，算符的另一边就与顶点位置为v0和v0+1
            lastspinop[s1] = v0 + 2
            lastspinop[s2] = v0 + 3
        else
            # 当前没有算符，顶点列表赋值为 0
            vertexlist[v0:v0+3] = [0, 0, 0, 0]
        end
    end

    # ============================================================
    # 创建“时间”边界上的最后链接
    # ============================================================
    @inbounds for s1 = 1:nn # 周期边界，头尾便利所有自旋nn链接，将第一次出现与最后一次出现链接
        v1 = firstspinop[s1]
        if v1 != -1
            v2 = lastspinop[s1]
            vertexlist[v2] = v1
            vertexlist[v1] = v2
        end
    end

    # ============================================================
    # 扫描更新：环路更新过程（优化后）
    # ============================================================
    loops_done::Int64 = 0 # 初始化已完成的循环次数
    while loops_done < nn
        v0 = rand(1:4*mm)
        if vertexlist[v0] < 1  # 跳过空算符
            continue
        end
        v1 = v0
        while true
            @inbounds begin
                # 使用局部变量缓存重复计算
                local_leg = v1 % 4
                enter_leg = (local_leg == 0 ? 4 : local_leg)
                nei_leg = (enter_leg in (1, 3)) ? (enter_leg + 1) : (enter_leg - 1)
                # 用 (v1+3)÷4 替代 cld(v1,4)
                bond_n = (v1 + 3) ÷ 4
                op = opstring[bond_n]
                b = fld(op, 2)
                # 预先缓存常用值
                local_spin1 = spin_mm[enter_leg, bond_n]
                local_spin2 = spin_mm[nei_leg, bond_n]
                site_enter = (enter_leg % 2 == 0 ? 2 : 1)
                site_nei = (nei_leg % 2 == 0 ? 2 : 1)
                local_offdig = J / 2
                local_dif = abs(Δ / 4) + Δ / 4 + W / 2
                local_upprob = abs(Δ / 4) - Δ / 4 + W / 2
                local_down = local_upprob
                local_prob = rand()

                if local_spin1 == local_spin2  # 自旋相同
                    if local_spin1 == 1
                        local_dif_h = local_dif - (h[bsites[site_enter, b]] - h[bsites[site_nei, b]]) / 4
                        local_upprob_h = local_upprob + (h[bsites[site_enter, b]] + h[bsites[site_nei, b]]) / 4
                        local_sum = local_dif_h + local_offdig + local_upprob_h
                        prob_back = local_upprob_h / local_sum
                        prob_through = local_dif_h / local_sum
                        if local_prob < prob_back
                            v1 = vertexlist[v1]
                        elseif local_prob < prob_back + prob_through
                            v2 = through(v1, spin_mm, enter_leg, bond_n)
                            v1 = vertexlist[v2]
                        else
                            v2 = crossthrough(v1, spin_mm, enter_leg, bond_n)
                            opstring[(v1+3)÷4] = xor(opstring[(v1+3)÷4], 1)
                            v1 = vertexlist[v2]
                        end
                    else
                        local_dif_h = local_dif + (h[bsites[site_enter, b]] - h[bsites[site_nei, b]]) / 4
                        local_down_h = local_down - (h[bsites[site_enter, b]] + h[bsites[site_nei, b]]) / 4
                        local_sum = local_dif_h + local_offdig + local_down_h
                        prob_back = local_down_h / local_sum
                        prob_through = local_dif_h / local_sum
                        if local_prob < prob_back
                            v1 = vertexlist[v1]
                        elseif local_prob < prob_back + prob_through
                            v2 = through(v1, spin_mm, enter_leg, bond_n)
                            v1 = vertexlist[v2]
                        else
                            v2 = crossthrough(v1, spin_mm, enter_leg, bond_n)
                            opstring[(v1+3)÷4] = xor(opstring[(v1+3)÷4], 1)
                            v1 = vertexlist[v2]
                        end
                    end
                elseif mod(op, 2) == 0  # 自旋相反且为对角算符
                    if local_spin1 == 1
                        local_dif_h = local_dif + (h[bsites[site_enter, b]] - h[bsites[site_nei, b]]) / 4
                        local_down_h = local_down - (h[bsites[site_enter, b]] + h[bsites[site_nei, b]]) / 4
                        local_sum = local_dif_h + local_offdig + local_down_h
                        prob_back = local_dif_h / local_sum
                        prob_through = local_down_h / local_sum
                        if local_prob < prob_back
                            v1 = vertexlist[v1]
                        elseif local_prob < prob_back + prob_through
                            v2 = through(v1, spin_mm, enter_leg, bond_n)
                            v1 = vertexlist[v2]
                        else
                            v2 = crossback(v1, spin_mm, enter_leg, bond_n)
                            opstring[(v1+3)÷4] = xor(opstring[(v1+3)÷4], 1)
                            v1 = vertexlist[v2]
                        end
                    else
                        local_dif_h = local_dif - (h[bsites[site_enter, b]] - h[bsites[site_nei, b]]) / 4
                        local_upprob_h = local_upprob + (h[bsites[site_enter, b]] + h[bsites[site_nei, b]]) / 4
                        local_sum = local_dif_h + local_offdig + local_upprob_h
                        prob_back = local_dif_h / local_sum
                        prob_through = local_upprob_h / local_sum
                        if local_prob < prob_back
                            v1 = vertexlist[v1]
                        elseif local_prob < prob_back + prob_through
                            v2 = through(v1, spin_mm, enter_leg, bond_n)
                            v1 = vertexlist[v2]
                        else
                            v2 = crossback(v1, spin_mm, enter_leg, bond_n)
                            opstring[(v1+3)÷4] = xor(opstring[(v1+3)÷4], 1)
                            v1 = vertexlist[v2]
                        end
                    end
                else  # 非对角算符情况
                    if local_spin1 == 1
                        local_dif_h = local_dif - (h[bsites[site_enter, b]] - h[bsites[site_nei, b]]) / 4
                        local_down_h = local_down - (h[bsites[site_enter, b]] + h[bsites[site_nei, b]]) / 4
                        local_sum = local_dif_h + local_offdig + local_down_h
                        prob_back = local_offdig / local_sum
                        prob_crossback = local_dif_h / local_sum
                        if local_prob < prob_back
                            v1 = vertexlist[v1]
                        elseif local_prob < prob_back + prob_crossback
                            v2 = crossback(v1, spin_mm, enter_leg, bond_n)
                            opstring[(v1+3)÷4] = xor(opstring[(v1+3)÷4], 1)
                            v1 = vertexlist[v2]
                        else
                            v2 = crossthrough(v1, spin_mm, enter_leg, bond_n)
                            opstring[(v1+3)÷4] = xor(opstring[(v1+3)÷4], 1)
                            v1 = vertexlist[v2]
                        end
                    else
                        local_dif_h = local_dif + (h[bsites[site_enter, b]] - h[bsites[site_nei, b]]) / 4
                        local_upprob_h = local_upprob + (h[bsites[site_enter, b]] + h[bsites[site_nei, b]]) / 4
                        local_sum = local_dif_h + local_offdig + local_upprob_h
                        prob_back = local_offdig / local_sum
                        prob_crossback = local_dif_h / local_sum
                        if local_prob < prob_back
                            v1 = vertexlist[v1]
                        elseif local_prob < prob_back + prob_crossback
                            v2 = crossback(v1, spin_mm, enter_leg, bond_n)
                            opstring[(v1+3)÷4] = xor(opstring[(v1+3)÷4], 1)
                            v1 = vertexlist[v2]
                        else
                            v2 = crossthrough(v1, spin_mm, enter_leg, bond_n)
                            opstring[(v1+3)÷4] = xor(opstring[(v1+3)÷4], 1)
                            v1 = vertexlist[v2]
                        end
                    end
                end
            end
            if v1 == v0  # 回到起点则结束当前环
                break
            end
        end

        loops_done += 1  # 完成一个环的更新
    end

    # ============================================================
    # 自旋更新：根据顶点标记更新自旋状态，以实现环翻转
    # firstspinop[i] 是自旋 ( s_i ) 在顶点链表中的第一个顶点编号，是自旋与顶点链表的唯一入口点。因此只需要考虑链表开始的地方的自旋
    # ============================================================
    @inbounds for i = 1:nn
        if firstspinop[i] != -1 # 这个自旋被加入过链表，第一次出现对应的顶点编号，这里不可能出现3，4号这种顶点
            bond_n = cld(firstspinop[i], 4) # 该顶点对应的bond编号
            enter_leg = mod(firstspinop[i], 4) == 0 ? 4 : mod(firstspinop[i], 4) # 计算进入的腿编号，1-4，mod取模

            spin[i] = spin_mm[enter_leg, bond_n] # 更新自旋状态
        else
            # 对于没有算符作用的孤立自旋，以 50% 概率翻转其自旋，这个自旋在m长度的算符串中没有任何算符作用，因此可以随意翻转
            if rand() < 0.5
                spin[i] = -spin[i]
            end
        end
    end

    # println(spin)
end  #end function loopupdate

function through(v1::Int64, spin_mm::Matrix{Int64}, enter_leg::Int64, bond_n::Int64)
    spin_mm[enter_leg, bond_n] = -spin_mm[enter_leg, bond_n] # 当前顶点作用自旋反向

    v2 = enter_leg in (1, 2) ? v1 + 2 : v1 - 2 # 目标位点
    out_leg = mod(v2, 4) == 0 ? 4 : mod(v2, 4) # 目标位点的腿编号，1-4，mod取模

    spin_mm[out_leg, bond_n] = -spin_mm[out_leg, bond_n] # 目标顶点作用自旋反向
    return v2
end

function crossback(v1::Int64, spin_mm::Matrix{Int64}, enter_leg::Int64, bond_n::Int64)
    spin_mm[enter_leg, bond_n] = -spin_mm[enter_leg, bond_n] # 当前顶点作用自旋反向

    v2 = enter_leg in (1, 3) ? v1 + 1 : v1 - 1 # 目标位点
    out_leg = mod(v2, 4) == 0 ? 4 : mod(v2, 4) # 目标位点的腿编号，1-4，mod取模

    spin_mm[out_leg, bond_n] = -spin_mm[out_leg, bond_n] # 目标顶点作用自旋反向
    return v2
end

function crossthrough(v1::Int64, spin_mm::Matrix{Int64}, enter_leg::Int64, bond_n::Int64)
    spin_mm[enter_leg, bond_n] = -spin_mm[enter_leg, bond_n] # 当前顶点作用自旋反向

    v2 = 8 * (bond_n - 1) + 5 - v1 # 目标位点
    out_leg = mod(v2, 4) == 0 ? 4 : mod(v2, 4) # 目标位点的腿编号，1-4，mod取模

    spin_mm[out_leg, bond_n] = -spin_mm[out_leg, bond_n] # 目标顶点作用自旋反向
    return v2
end

function adjustcutoff(nh::Int64, mm::Int64, step::Int64,
    opstring::Vector{Int64},
    vertexlist::Vector{Int64},
    spin_mm::Matrix{Int64})  # 4 × mm 的矩阵

    mmnew::Int64 = nh + fld(nh, 3)
    if mmnew > mm
        resize!(opstring, mmnew)
        opstring[mm+1:mmnew] .= 0
        resize!(vertexlist, 4 * mmnew)

        # 手动扩展二维数组（4 × mm 的结构）
        spin_mm_new = zeros(Int64, 4, mmnew)         # 新的是 4 × mmnew
        spin_mm_new[:, 1:mm] .= spin_mm              # 拷贝原数据
        spin_mm = spin_mm_new                        # 替换
        mm = mmnew

        open("XXZ_2D_afm_random_field/cut.text", "a") do f
            println(f, "step: ", step, " Cut-off L: ", mmnew)
        end
    end

    return mm, spin_mm
end

end # module SSEUpdate
