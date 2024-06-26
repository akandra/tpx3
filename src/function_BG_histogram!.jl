include("find_TOF_interval!.jl")

function BG_histogram!(BGhists,
                     p::pars, 
                     out1::Vector{Float32}, 
                     out2::Vector{Float32}, 
                     out3::Vector{Float32}, 
                     out4::Vector{UInt16}, 
                     out5::Vector{UInt8}, 
                     out6::Vector{UInt8},
                     n_out::Int, 
                     ROI_sig::Matrix{Int16},
                     ROI_bkg::Matrix{Int16})

    if isnothing(ROI_sig) || isnothing(ROI_bkg)

        println("background histogram not constructed: no ROI defined")

    else

        ROI_for_hist = [[Int16(0) Int16(257)] ; [ROI_bkg[2,1] ROI_sig[2,2]]] # syntax: [[xmin, xmax]; [ymin, ymax]]
        
        #println("Getting ROI data...")
        idx_roi      = get_ROIdata(out1, out2, out3, out4, out5, out6, n_out, ROI_for_hist)
        
        #println("Constructing BGhists array...")
        for i in idx_roi
            if out1[i] > 0 && out1[i] < p.first_seconds
                
                which_hist = find_TOF_interval!(out3[i], p)
                
                if which_hist > -1
                    BGhists[which_hist][Int(out5[i])] += 1
                end
                
            end
        end

    end
    
end