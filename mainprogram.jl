module SSEMain
export mainprogram

include("lattice.jl")
include("update.jl")
include("measure.jl")

using .SSELattice, .SSEUpdate, .SSEMeasure

# 清空或创建测量结果文件
open("XXZ_2D_afm_random_field/res_2D.dat", "w") do f
    # 创建并清空文件
end

# 清空或创建测量结果文件
open("XXZ_2D_afm_random_field/cut.text", "w") do f
    # 创建并清空文件
end

function mainprogram(lx::Int64, ly::Int64, beta::Float64, nbins::Int64, msteps::Int64, isteps::Int64)
    # lx: 链长, beta: 逆温, nbins: bin 数量, msteps: 测量步数, isteps: 平衡步数
    nn::Int64 = lx * ly         # 1D 链，总自旋数即链长
    nb::Int64 = 2 * nn        # 键总数（每个自旋参与两个键）
    mm::Int64 = max(4, fld(nn, 4))  # 算符串截断长度，fld表示 向下取整的整数除法
    nh::Int64 = 0             # 当前对角算符数
    Δ::Float64 = -1.0          # z方向自旋耦合常数
    J::Float64 = 1.0          # x，y方向自旋耦合常数
    W::Float64 = 0.0          # 无序
    # h = W * (2 * rand(lx) .- 1) # 随机场强度，rand(lx)
    h = fill(W, nn) # 恒定外场
    # h = [0.0032726873244453003, 0.04815983495745333, -0.03472934390892135, -0.106279585602918, 0.150402539024184, -0.14377758767769763, 0.13826873090363528, -0.19771812466580052]

    # 初始化自旋状态，取值为 -1 或 1（通过 rand(-1:2:1) 随机取值），生成一个长度为 lx 的向量，元素是服从**[0, 1)** 均匀分布的随机数。2 * rand(lx) .- 1将上一步的向量映射到区间 [-1, 1)：
    spin = Array{Int64}(undef, nn)
    @inbounds for i = 1:nn
        spin[i] = rand(-1:2:1)
    end

    spin_mm = zeros(Int64, 4, mm)  # 行是 mm，列是4
    # 记录每个键的算符种类
    bsites = Array{Int64}(undef, 2, nb)
    opstring = zeros(Int64, mm)
    vertexlist = Array{Int64}(undef, 4 * mm) # 记录链接表，记录所有顶点的连接情况

    bsites = SSELattice.makelattice(nn, lx, ly, bsites)  # 调用 SSELattice 模块中的函数

    # 平衡过程，进行 isteps 步对角和循环更新
    for i = 1:isteps
        nh = SSEUpdate.diagonalupdate(nh, mm, opstring, spin, spin_mm, bsites, nb, beta, Δ, J, h, W)
        SSEUpdate.loopupdate(mm, nn, opstring, bsites, vertexlist, spin, spin_mm, Δ, J, h, W)
        # SSEUpdate.loopupdate(mm, nn, opstring, bsites, vertexlist, spin)
        mm, spin_mm = SSEUpdate.adjustcutoff(nh, mm, i, opstring, vertexlist, spin_mm)
    end

    # 主测量过程，每个 bin 内运行 msteps 步更新和测量
    for j = 1:nbins
        enrg1::Float64 = 0.0
        enrg2::Float64 = 0.0
        amag::Float64 = 0.0
        amag2::Float64 = 0.0
        MRC::Float64 = 0.0
        for i = 1:msteps
            nh = SSEUpdate.diagonalupdate(nh, mm, opstring, spin, spin_mm, bsites, nb, beta, Δ, J, h, W,)
            SSEUpdate.loopupdate(mm, nn, opstring, bsites, vertexlist, spin, spin_mm, Δ, J, h, W)
            # SSEUpdate.loopupdate(mm, nn, opstring, bsites, vertexlist, spin)
            (enrg1, enrg2, amag, amag2, MRC) = SSEMeasure.measureobservables(mm, nh, nn, lx, ly, spin, bsites, opstring, enrg1, enrg2, amag, amag2, MRC)
        end
        SSEMeasure.writeresults(Δ, W, nn, nb, msteps, beta, enrg1, enrg2, amag, amag2, MRC)

        println(spin)
    end
end

# 参数读取部分，从 "read.in" 文件中读取模拟参数
if isfile("XXZ_2D_afm_random_field/read_2D.in")
    f = open("XXZ_2D_afm_random_field/read_2D.in", "r")
    lx = parse(Int64, split(readline(f))[1])
    ly = parse(Int64, split(readline(f))[1])
    beta = parse(Float64, split(readline(f))[1])
    nbins = parse(Int64, split(readline(f))[1])
    isteps = parse(Int64, split(readline(f))[1])
    msteps = parse(Int64, split(readline(f))[1])
    close(f)
else
    println("Warning: 'read.in' not found. Using default parameters.")
    lx, beta, nbins, isteps, msteps = 16, 1.0, 100, 10, 1000
end

@time mainprogram(lx, ly, beta, nbins, msteps, isteps)

end # module SSEMain