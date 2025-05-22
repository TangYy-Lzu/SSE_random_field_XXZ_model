module SSEMeasure
export measureobservables, writeresults

using LinearAlgebra        # 使用线性代数库

# 物理量测量函数和结果写入函数
# 增加详细注释帮助理解测量过程，包括初始自旋构型、算符作用以及最终结果的处理

function measureobservables(mm::Int64, nh::Int64, nn::Int64, lx::Int64, spin::Vector{Int64}, bsites::Matrix{Int64}, opstring::Vector{Int64}, enrg1::Float64, enrg2::Float64, amag::Float64, amag2::Float64, MRC::Float64)
    # 计算系统自旋配置下的物理量
    # mm: 当前算符字符串长度, nh: 对角算符数目, nn: 自旋数, lx: 格子宽度
    # spin: 自旋数组, bsites: 键与自旋的对应, opstring: 算符字符串
    # enrg1, enrg2, amag2, ususc: 分别表示能量、能量平方、磁化平方和磁化易感度累加量

    am::Int64 = 0
    # 初始阶段计算交错自旋的磁化（staggered magnetization）
    for i = 1:nn
        am += spin[i] * (-1)^(i - 1)
    end
    am = fld(am, 2)  # 初始自旋总和除以2

    am_mean::Float64 = 0.0  # 平均交错磁化
    am2::Float64 = 0.0  # 磁化平方，因为第一行自旋与最后一行相同（周期边界）因此这一行无需考虑
    MRC_step::Float64 = 0.0  # 最大特征值

    M = zeros(Float64, nn, nn)
    # 遍历算符字符串，处理每个操作对磁化的影响，对所有虚时遍历
    for i = 1:mm
        op = opstring[i]
        if op < 0 # 判断为外场算符
            continue # 跳过外场算符
        elseif mod(op, 2) == 1 # 判断为非对角算符（off-diagonal operator）
            b = fld(op, 2)
            s1 = bsites[1, b]
            s2 = bsites[2, b]
            spin[s1] = -spin[s1]
            spin[s2] = -spin[s2]
            # 对于非对角操作，翻转相关自旋后调整交错磁化
            am = am + 2 * spin[s1] * (-1)^(mod(s1 - 1, lx))
            am_mean = am_mean + am
        end
        for i = 1:nn
            for j in i:nn
                M[i, j] += spin[i] * spin[j]
                if i != j
                    M[j, i] += M[i, j]
                end
            end
        end
        am2 = am2 + am^2  # 累加当前的磁化平方
    end
    # println(M / (4 * mm))
    MRC_step += eigmax(Symmetric(M)) / (4 * nn * mm) # 累加MRC
    # println("MRC_step: ", MRC_step)

    am_mean = am_mean / mm
    am2 = am2 / mm  # 求平均值
    enrg1 = enrg1 + nh
    enrg2 = enrg2 + nh^2
    MRC_step = MRC_step

    amag = amag + am_mean
    amag2 = amag2 + am2
    MRC = MRC + MRC_step  # 平均最大特征值
    # ususc = ususc + (sum(spin) / 2)^2  # 自旋总和平方贡献自磁化易感度
    return enrg1, enrg2, amag, amag2, MRC
end

function writeresults(Δ::Float64, W::Float64, nn::Int64, nb::Int64, msteps::Int64, beta::Float64, enrg1::Float64, enrg2::Float64, amag::Float64, amag2::Float64, MRC::Float64)
    # 对各物理量做步数归一化，然后写入结果文件
    enrg1 = enrg1 / msteps
    enrg2 = enrg2 / msteps
    amag = amag / msteps
    amag2 = amag2 / msteps
    MRC = MRC / msteps
    # ususc = ususc / msteps

    # 调整能量数据，将测量结果归一化为每个自旋的贡献
    enrg2 = (enrg2 - enrg1 * (enrg1 + 1.0)) / nn
    enrg1 = -(enrg1 / (beta * nn) - (abs(0.25 * Δ) + 0.5 * W) * nb / nn) # 因为我们给哈密顿量加上了一个0.25的因子以让对角项作用为0
    amag = 3.0 * amag / nn
    amag2 = 3.0 * amag2 / (nn^2) # 乘以三是因为三个方向的总贡献
    # ususc = beta * ususc / nn

    # 将归一化后的测量结果追加写入结果文件
    f = open("SSE-Heisenberg_1D/res_1D.dat", "a")
    println(f, enrg1, "   ", enrg2, "   ", amag, "   ", amag2, "   ", MRC)
    close(f)
    # 写入结果后，可在此处清零累加变量以便后续测量
end

end # module SSEMeasure
