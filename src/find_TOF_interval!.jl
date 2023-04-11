function find_TOF_interval!(tof_value::Float32, p::pars)
    # returns the index of the respective TOF interval if
    # lower border of TOF gate < tof_value <= upper border of TOF gate
    # else: returns nothing

    idx = nothing

    if p.tof_gates[1] <= tof_value <= p.tof_gates[end]
        i = 1
        while tof_value > p.tof_gates[i]
            i += 1
        end
        if mod(i, 2) == 0
            idx = Int(i/2)
        end
    end
    
    return(idx)
end