# lattice.jl

module SSELattice
export makelattice

# 构造一维链的函数
# 详细注释：解释如何根据 L 计算每个自旋的位置及生成对应的键，键看为哈密顿量作用的地方

function makelattice(nn::Int64, lx::Int64, ly::Int64, bsites::Matrix{Int64})
    @inbounds for y1 = 0:ly-1 # @inbounds不检查数组索引是否越界，提高性能
        @inbounds for x1 = 0:lx-1
            s::Int64 = 1 + x1 + y1 * lx  # 当前格点编号（从1开始）
            # 水平键：右侧的邻居，周期性边界条件
            x2 = mod(x1 + 1, lx)
            y2 = y1
            bsites[1, s] = s # 哈密顿量作用在的水平方向的两个格点的位置
            bsites[2, s] = 1 + x2 + y2 * lx

            # 垂直键：下方的邻居，周期性边界条件
            x2 = x1
            y2 = mod(y1 + 1, ly)
            bsites[1, s+nn] = s # 哈密顿量作用在的垂直方向的两个格点的位置，nn为格点数，前s个表示了水平位置的bond
            bsites[2, s+nn] = 1 + x2 + y2 * lx
        end
    end
    return bsites
end

end # module SSELattice