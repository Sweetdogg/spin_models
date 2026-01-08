using Statistics
using DelimitedFiles
using Random
using LinearAlgebra

mutable struct PottsSystem
    L::Int
    N::Int
    q::Int
    T::Float64
    spins::Vector{Int8}
    nbor::Matrix{Int32}
    p_add::Float64
end

function init_system(L::Int, T::Float64, q::Int=3)
    N = L * L
    spins = rand(Int8(0):Int8(q-1), N)
    nbor = zeros(Int32, 4, N)
    
    for i in 1:N
        x = (i - 1) % L
        y = (i - 1) ÷ L
        
        right = ((x + 1) % L) + y * L + 1
        left  = ((x - 1 + L) % L) + y * L + 1
        up    = x + ((y + 1) % L) * L + 1
        down  = x + ((y - 1 + L) % L) * L + 1
        
        nbor[1, i] = Int32(right)
        nbor[2, i] = Int32(up)
        nbor[3, i] = Int32(left)
        nbor[4, i] = Int32(down)
    end

    p_add = 1.0 - exp(-1.0 / T)
    return PottsSystem(L, N, q, T, spins, nbor, p_add)
end

function wolff!(sys::PottsSystem)
    N = sys.N
    q = sys.q
    spins = sys.spins
    nbor = sys.nbor
    p_add = sys.p_add
    
    seed_idx = rand(1:N)
    old_spin = spins[seed_idx]
    new_spin = Int8(mod(old_spin + rand(1:(q-1)), q))
    
    stack = Int32[seed_idx]
    spins[seed_idx] = new_spin
    
    @inbounds while !isempty(stack)
        current_idx = pop!(stack)
        
        for k in 1:4
            neighbor_idx = nbor[k, current_idx]
            if spins[neighbor_idx] == old_spin
                if rand() < p_add
                    spins[neighbor_idx] = new_spin
                    push!(stack, neighbor_idx)
                end
            end
        end
    end
end

# --- NEW: Energy Calculation ---
function measure_energy(sys::PottsSystem)
    E = 0.0
    N = sys.N
    spins = sys.spins
    nbor = sys.nbor
    
    @inbounds for i in 1:N
        right = nbor[1, i]
        up    = nbor[2, i]
        
        if spins[i] == spins[right]
            E -= 1.0
        end
        if spins[i] == spins[up]
            E -= 1.0
        end
    end
    return E
end

# --- UPDATED: Measure everything ---
function measure_observables(sys::PottsSystem)
    # 1. Magnetization
    mx = 0.0
    my = 0.0
    factor = 2.0 * π / sys.q
    
    @inbounds for s in sys.spins
        angle = s * factor
        s_sin, s_cos = sincos(angle)
        mx += s_cos
        my += s_sin
    end
    
    m2 = (mx^2 + my^2) / (sys.N^2)
    m_abs = sqrt(m2)
    
    # 2. Energy
    E = measure_energy(sys)
    
    return E, m_abs, m2
end

function run_simulation(L, T, mc_eq, bins, mc_av)
    sys = init_system(L, T, 3)
    N = sys.N
    
    # Equilibration
    for _ in 1:mc_eq
        wolff!(sys)
    end
    
    # Arrays to store results per bin
    res_binder = zeros(Float64, bins)
    res_chi    = zeros(Float64, bins)
    res_cv     = zeros(Float64, bins)
    res_m2     = zeros(Float64, bins)
    res_E      = zeros(Float64, bins)
    
    # Sampling
    for b in 1:bins
        E_sum = 0.0
        E2_sum = 0.0
        m_sum = 0.0
        m2_sum = 0.0
        m4_sum = 0.0
        
        for _ in 1:mc_av
            wolff!(sys)
            
            E, m, m2 = measure_observables(sys)
            
            E_sum += E
            E2_sum += E^2
            m_sum += m
            m2_sum += m2
            m4_sum += m2^2
        end
        
        # Averages for this bin
        avg_E = E_sum / mc_av
        avg_E2 = E2_sum / mc_av
        avg_m = m_sum / mc_av
        avg_m2 = m2_sum / mc_av
        avg_m4 = m4_sum / mc_av
        
        
        
        # 1. Binder Cumulant
        res_binder[b] = avg_m4 / (avg_m2^2)
        
        # 2. Susceptibility: chi = N/T * (<m^2> - <|m|>^2)
        res_chi[b] = (N / T) * (avg_m2 - avg_m^2)
        
        # 3. Heat Capacity: Cv = ( <E^2> - <E>^2 ) / (N * T^2)
        res_cv[b] = (avg_E2 - avg_E^2) / (N * T^2)
        
        # 4. Magnetization squared
        res_m2[b] = avg_m2
        
        # 5. Energy per site
        res_E[b] = avg_E / N
    end
    
    return res_binder, res_chi, res_cv, res_m2, res_E
end

# --- Main Execution ---

if isfile("read.in")
    params = readdlm("read.in")
    L_in = Int(params[1])
    T_in = Float64(params[2])
    mc_eq_in = Int(params[3])
    bins_in = Int(params[4])
    mc_av_in = Int(params[5])
else
    L_in = 64
    T_in = 1.0
    mc_eq_in = 2000
    bins_in = 10
    mc_av_in = 5000
end

binders, chis, cvs, m2s, es = run_simulation(L_in, T_in, mc_eq_in, bins_in, mc_av_in)

function get_stats(data)
    return mean(data), std(data) / sqrt(length(data))
end

u_mean, u_err = get_stats(binders)
chi_mean, chi_err = get_stats(chis)
cv_mean, cv_err = get_stats(cvs)
m2_mean, m2_err = get_stats(m2s)
e_mean, e_err = get_stats(es)

# Output raw bin data (Optional, maybe just output averages)
open("out.dat", "a") do io
    for i in 1:bins_in
        # Format: L T Bin Binder Chi Cv M2 E
        println(io, "$L_in $T_in $i $(binders[i]) $(chis[i]) $(cvs[i]) $(m2s[i]) $(es[i])")
    end
end

# Output averaged result
open("out_avg.dat", "w") do io
    # Format: L T U U_err Chi Chi_err Cv Cv_err M2 M2_err E E_err
    println(io, "$L_in $T_in $u_mean $u_err $chi_mean $chi_err $cv_mean $cv_err $m2_mean $m2_err $e_mean $e_err")
end

# flush输出文件？
