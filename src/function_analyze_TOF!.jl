function analyze_TOF!(  tof_hist::Vector{Int32}, 
                        p::pars, 
                        beam_time::Vector{Float32}, 
                        bl_delay::Vector{Float32}, 
                        tof::Vector{Float32}, 
                        tot::Vector{UInt16}, 
                        x::Vector{UInt8}, 
                        y::Vector{UInt8}, 
                        n_out::Int, 
                        ROI::Matrix{Int16})

    idx_list = get_ROIdata(beam_time, bl_delay, tof, tot, x, y, n_out, ROI)
    
    for i in idx_list
        idx = Int(floor(tof[i] / p.tof_bin) + 1)  # idx 1 corresponds to bin with time=0 to TOFbin
        if idx > 1 && idx <= length(tof_hist)  # ignore scattered light ions in bin 1
            tof_hist[idx] += 1
        end
    end
    
end