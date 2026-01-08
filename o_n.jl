using Statistics
using DelimitedFiles

struct SimParams
    l::Int
    temp::Float64
    mc_eq::Int
    bins::Int
    mc_av::Int
    on::Int
end

function metropolis(spn, nbor, temp)
    on, n = size(spn)
    enbor = zeros(Float64, on)
    for i in 1:n
        vt = randn(on)
        vt /= sqrt(sum(vt.^2))
        idx = rand(1:n)
        s = spn[:, idx]
        enbor = sum(spn[:, nbor[:, idx]], dims=2)
        ei =  (enbor's)[1]
        ef =  (enbor'vt)[1]
        pj = exp((ef - ei) / temp) # pj怎么可能是个矩阵
        if rand() < pj
            spn[:, idx] = vt
        end
    end
    return spn
end

function wolff(spn, nbor, temp)
    on, n = size(spn)
    idx = rand(1:n)
    s = spn[:, idx]
    r = randn(on)
    r /= sqrt(sum(r.^2))
    c = s'r
    cluster = falses(n)
    cluster[idx] = true
    stack = [idx]
    while !isempty(stack)
        current = pop!(stack)
        for j in nbor[:, current]
            if !cluster[j]
                sj = spn[:, j]
                pj = 1 - exp((2 * (sj'r) * (s'r)) / temp)
                if rand() < pj
                    cluster[j] = true
                    push!(stack, j)
                end
            end
        end
    end
    for i in 1:n
        if cluster[i]
            si = spn[:, i]
            spn[:, i] -= 2 * (si'r) * r
        end
    end
    return spn
end

function run(isteps, bins, steps, spn, temp, nbor)
    mag2 :: Float64 = 0
    mag4 :: Float64 = 0
    binder = zeros(Float64, bins)

    for _ in 1:isteps
        spn = wolff(spn, nbor, temp)
        # spn = metropolis(spn, nbor, temp)
    end
    for i in 1:bins
        mag2 = 0; mag4 = 0
        for _ in 1:steps
            # spn = metropolis(spn, nbor, temp)
            spn = wolff(spn, nbor, temp)
            m_vec = sum(spn, dims=2)
            m2 = sum(m_vec.^2)
            mag2 += m2
            mag4 += m2^2
        end
        binder[i] = mag4/mag2^2 * steps
    end
    return binder
end

l, temp, mc_eq, bins, mc_av, on = readdlm("read.in")
pa = SimParams(Int(l), Float64(temp), Int(mc_eq), Int(bins), Int(mc_av), Int(on))

n = pa.l^2
spn = randn(pa.on, n)
spn ./= sqrt.(sum(spn.^2, dims=1))
nbor = zeros(Int32, 4, n)

for i in 1:n
    x = (i - 1) % pa.l
    y = (i - 1) ÷ pa.l
    nbor[1, i] = ((x + 1) % pa.l) + y * pa.l + 1
    nbor[2, i] = ((x - 1 + pa.l) % pa.l) + y * pa.l + 1
    nbor[3, i] = x + ((y + 1) % pa.l) * pa.l + 1
    nbor[4, i] = x + ((y - 1 + pa.l) % pa.l) * pa.l + 1
end
binders = run(pa.mc_eq, pa.bins, pa.mc_av, spn, pa.temp, nbor)
binder = mean(binders)
error_val = std(binders) / sqrt(pa.bins)
# 原始数据输出
io1 = open("out.dat", "a")
for i in 1:pa.bins
    println(io1, "$(pa.l) $(pa.temp) $(binders[i])")
end
close(io1)

# 处理数据输出
io2 = open("res.dat", "a")
println(io2, "$(pa.l) $(pa.temp) $(binder) $(error_val)")
close(io2)

# 能否创建一个结构体，包括spn, nbor, on, n等参数以及蒙卡更新的方法
# 注意！的用法
# 矩阵的维度问题