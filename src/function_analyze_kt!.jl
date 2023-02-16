function analyze_kt!(kt_hist::Vector{Vector{Int32}}, 
                     p::pars, 
                     beam_time::Vector{Float32}, 
                     bl_delay::Vector{Float32}, 
                     tof::Vector{Float32}, 
                     tot::Vector{UInt16}, 
                     x::Vector{UInt8}, 
                     y::Vector{UInt8}, 
                     n_out::Int, 
                     ROI::Matrix{Int16})

    idx_roi = get_ROIdata(beam_time, bl_delay, tof, tot, x, y, n_out, ROI)

    dummy_counter = 0
    
    for i in idx_roi

        if beam_time[i] > p.kt_t0
            if tof[i] > p.tof_gate_min && tof[i] < p.tof_gate_max
                idx_kt = Int(floor( (beam_time[i] - p.kt_t0) / p.kt_length + 1 ))
                
                if idx_kt > length(kt_hist)
                    #println("\t \t \t hist switch occurs at index ", i, " and beam time ", beam_time[i])
                    push!(kt_hist, zeros(Int32, p.kt_nbins))
                end
                
                idx = Int(floor(bl_delay[i] / p.kt_bin) + 1)     # idx 1 corresponds to bin with time=0 to TOFbin
                
                if idx > 0 && idx <= length(kt_hist[idx_kt])
                    kt_hist[idx_kt][idx] += 1
                end
                
            end
        end
    end

    
end