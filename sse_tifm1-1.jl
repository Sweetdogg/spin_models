# 尺寸为 L 的一维横场伊辛模型的 SSE 模拟
# 参考文献：A. W. Sandvik, ”Stochastic series expansion method for quantum Ising models with arbitrary interactions”, Phys. Rev. E 68, 056701 (2003).
mutable struct M
    L::Int
    h::Float64  # 横场强度
    beta::Float64 
    m::Int # 序列长度
    n::Int # 棒数
    op::Array{Int,2} # 操作符序列，2行m列，第一行表示类型，第二行表示作用的自旋或键 
    s::Array{Int,1} # 自旋态
    b::Array{Int,2} # 键表，2行L列
    v::Array{Int,1} # 链表
    dat::Array{Float64,1}
end

init(beta, L, h) = M(L, h, beta, L^2, 0, zeros(Int,2,L^2), rand([-1,1],L), [1:L [2:L; 1]]', zeros(Int,4*L^2), zeros(7))

function run!(sys::M, meas=false, cutoff=true)
    N, nb, beta = sys.L, sys.L, sys.beta
    h = sys.h
    nid = sys.m - sys.n
    pp = (h + 2) * N * beta
    # --- Diagonal Update ---
    for i in 1:sys.m
        o, r = sys.op[:,i], rand()
        if o[1] == 0 # 插入
            if r * nid < pp
                id = rand(1:nb) 
                if rand() < h / (h + 2) # 横场
                    # sys.op[:,i] = [-1; id]
                    sys.op[1,i], sys.op[2,i] = -1, id
                    sys.n += 1; nid -= 1
                else
                    s1, s2 = sys.b[1,id], sys.b[2,id]
                    if sys.s[s1] == sys.s[s2] # 铁磁
                        # sys.op[:,i] = [s1; s2]
                        sys.op[1,i], sys.op[2,i] = s1, s2
                        sys.n += 1; nid -= 1
                    end
                end
            end
        elseif o[1] == -1 || o[1] > 0 # 移除
            if r < (nid + 1) / pp
                sys.op[:,i] = [0; 0]; sys.n -= 1; nid += 1
            end
        else # 翻转
            sys.s[o[2]] *= -1
        end
    end

    # 构建链表
    frst, last = fill(-1, N), fill(-1, N)
    vrtx = sys.v; vrtx .= 0
    for (i,o) in enumerate(eachcol(sys.op))
        o[1] == 0 && continue
        v0 = 4*(i-1) + 1
        if o[1] < 0 # transverse field
            s1 = o[2]
            v1 = last[s1]
            if frst[s1] != -1
                vrtx[v0], vrtx[v1] = v1, v0
            else
                frst[s1] = v0
            end
            last[s1] = v0+2          
        else # bond
            s1, s2 = o
            v1, v2 = last[s1], last[s2]
            if v1 != -1
                vrtx[v0], vrtx[v1] = v1, v0
            else
                frst[s1] = v0
            end
            if v2 != -1
                vrtx[v0+1], vrtx[v2] = v2, v0+1
            else
                frst[s2] = v0+1
            end
            last[s1], last[s2] = v0+2, v0+3
        end
    end
    for s in 1:N # 闭合链表
        if frst[s] != -1
            vrtx[last[s]], vrtx[frst[s]] = frst[s], last[s]
        end
    end

    # for i in 1:S.m
    #     println(i, ", ", S.op[:,i], ", ", S.v[4*(i-1)+1:4*i])
    # end

    # Classical cluster update 不改变横场项 
    # quantum cluster update（待优化）
    for v0 in 1:2:4*sys.m
        vrtx[v0] < 1 && continue
        cluster = Set{Int}()
        # cluster = Int[] --- ERROR IGNORE ---
        stack = [v0]
        
        # 构造集群
        while !isempty(stack)
            v1 = pop!(stack)
            # println("处理顶点 ", v1)
            v1 > 0 || continue
            op = sys.op[:, (v1+3) ÷ 4]
            push!(cluster, v1)
            # 遍历相连顶点
            if op[1] < 0 # transverse field
                v2 = vrtx[v1]
                push!(stack, v2)
                vrtx[v1] = 0
            else # bond
                # v1 = v1 ÷ 4
                v1 = (v1+3) ÷ 4
                v1 = 4*(v1-1) + 1
                for v2 in v1:v1+3
                    push!(stack, vrtx[v2])
                    vrtx[v2] = 0
                end
            end
        end
        
        # println(cluster)

        # 翻转集群
        if rand() < 0.5 
            for v in cluster
                vrtx[v] = -1
                s1 = (v+3) ÷ 4
                if sys.op[1, s1] < 0
                    sys.op[1, s1] = xor(sys.op[1, s1], 1)
                end
            end
        end
    end

    # 更新自旋
    for i in 1:sys.L
        fi = frst[i]
        if fi != -1
            if vrtx[fi] == -1
                sys.s[i] *= -1
            end
        else
            rand() < 0.5 && (sys.s[i] *= -1)
        end
    end


    # 调整截断
    if cutoff
        new_m = sys.n + sys.n÷3
        if new_m > sys.m
            sys.op = hcat(sys.op, zeros(Int, 2, new_m - sys.m))
            sys.m = new_m
            resize!(sys.v, 4*sys.m)
        end
    end

    if meas
        sys.dat[1] += sys.n
        # nflip = count(o -> o[1] == -2, eachcol(sys.op))
        # sys.dat[2] += nflip
        mag2 = (sum(sys.s)/Float64(sys.L))^2
        sys.dat[3] += sqrt(mag2)
        sys.dat[4] += mag2
        sys.dat[5] += mag2^2
    end
    
    # return sys.v
end
#  Bootstrap 和 Jackknife 方法用于误差估计
L = 16
beta = L
bins = 10
eq = 10000
mcs = 10000
S = init(beta, L, 1.0)
foreach(_->run!(S), 1:eq)
for _ in 1:bins
    S.dat .= 0.0
    foreach(_->run!(S, true, false), 1:mcs)
    S.dat ./= mcs
    E = -S.dat[1] / S.beta + (1 + S.h) * S.L
    M2 = S.dat[4]
    R = S.dat[5] / (M2^2)
    chi = beta * L *(M2 - (S.dat[3])^2)
    @show L, beta, S.h, E / L, chi, R, M2
end