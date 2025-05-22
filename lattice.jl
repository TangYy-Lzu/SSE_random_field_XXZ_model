# lattice.jl

module SSELattice
export makelattice

# 构造一维链的函数
# 详细注释：解释如何根据 L 计算每个自旋的位置及生成对应的键，键看为哈密顿量作用的地方

function makelattice(nn::Int64, L::Int64, bsites::Matrix{Int64})
    for s = 1:nn
        bsites[1, s] = s
        bsites[2, s] = s == nn ? 1 : s + 1
    end
    return bsites
end

end # module SSELattice
