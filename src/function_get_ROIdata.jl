function get_ROIdata(   beam_time::Vector{Float32}, 
                        bl_delay::Vector{Float32}, 
                        tof::Vector{Float32}, 
                        tot::Vector{UInt16}, 
                        x::Vector{UInt8}, 
                        y::Vector{UInt8},
                        n_out::Int, 
                        ROI::Matrix{Int16})

    if isnothing(ROI)

        idx_list = 1:n_out
        stop # should never happen
    else

        idx_list = zeros(Int,n_out)
        idx = [0]
        for i in 1:n_out
            if x[i] > ROI[1,1] && x[i] < ROI[1,2] && y[i] > ROI[2,1] && y[i] < ROI[2,2]
                #push!(idx_list, i) 
                idx[1] += 1
                idx_list[idx[1]] = i
            end
        end
    
    end

    return resize!(idx_list,idx[1])
end
