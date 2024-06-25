include("find_TOF_interval!.jl")

function realtime_analysis!(realtime_trace, p, ROI, beam_time, bl_delay, tof, tot, x, y, n_out)
    
    idx_roi           =  get_ROIdata(beam_time, bl_delay, tof, tot, x, y, n_out, ROI)
    
    for i in idx_roi
        
        if beam_time[i] > p.realtime_t0 && beam_time[i] < p.realtime_tmax

            which_kt      = find_TOF_interval!(tof[i], p)


            if which_kt > -1

                lasertime = beam_time[i] + bl_delay[i]

                binnr     = Int(floor( (lasertime - p.realtime_t0) / p.realtime_bin) + 1)

                if binnr > 0 && binnr <= size(realtime_trace[which_kt])[1]
                    ion_velocity  = abs(x[i] - p.image_x_laser) / tof[i] / cos(3.14159 * p.measurement_angle / 180) / p.pixel_per_mm * 1e-3 # m/s
                    realtime_trace[which_kt][binnr] += 1 * Int64(round(ion_velocity, digits=0))
                end
                
            end
        end
    end
end