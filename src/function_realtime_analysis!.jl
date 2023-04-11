include("find_TOF_interval!.jl")

function realtime_analysis!(realtime_trace, p, ROI, beam_time, bl_delay, tof, tot, x, y, n_out)
    
    idx_roi           =  get_ROIdata(beam_time, bl_delay, tof, tot, x, y, n_out, ROI)
    
    for i in idx_roi
        
        if beam_time[i] > p.realtime_t0 && beam_time[i] < p.realtime_tmax

            which_kt      = find_TOF_interval!(tof[i], p)


            if typeof(which_kt) != Nothing

                lasertime = beam_time[i] + bl_delay[i]

                binnr     = Int(floor( (lasertime - p.realtime_t0) / p.realtime_bin) + 1)

                if binnr > 0 && binnr <= size(realtime_trace[which_kt])[1]
                    realtime_trace[which_kt][binnr] += 1
                end
                
            end
        end
    end
end