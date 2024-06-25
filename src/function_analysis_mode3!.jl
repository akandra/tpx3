include("find_TOF_interval!.jl")

function analysis_mode3!(tof_hist, kt_histogr, kt_histogr_first_sec, p, ROI, beam_time, bl_delay, tof, tot, x, y, n_out)
    
    idx_roi       =  get_ROIdata(beam_time, bl_delay, tof, tot, x, y, n_out, ROI)

    # NOW, LOOP ONCE OVER THE EVENTS IN THE ROI
    
    #Threads.@threads for i in idx_roi 
    # when using this, the number of predefined trace_chunks (default: 20) must be sufficient, i.e. no push! should be used!
        
    for i in idx_roi
        
        # 1. GENERATE TOF SPECTRUM
        if beam_time[i] > p.mode3_t0 && beam_time[i] < p.mode3_tmax
            tofidx   =   Int(floor(tof[i] / p.tof_bin) + 1)    
            if tofidx > 1  &&  tofidx <= length(tof_hist)
                tof_hist[tofidx] += 1
            end
        end
            
        # 2. GENERATE KINETIC TRACE FOR FIRST SECONDS
        if beam_time[i] > 0 && beam_time[i] < p.first_seconds
            
            which_kt = find_TOF_interval!(tof[i], p)
            
            if which_kt > -1
                binnr = Int(floor(bl_delay[i] / p.kt_bin) + 1)
                
                if binnr > 0 && binnr <= length(kt_histogr_first_sec[which_kt])
                    kt_histogr_first_sec[which_kt][binnr] += 1
                end
            end
        end
        
        # 3. GENERATE KINETIC TRACES IN CHUNKS FOR TIMES > mode3_t0
        if beam_time[i] > p.mode3_t0 && beam_time[i] < p.mode3_tmax
            
            which_kt = find_TOF_interval!(tof[i], p)
            
            if which_kt > -1
                kt_chunk = Int(floor( (beam_time[i] - p.mode3_t0) / p.kt_length + 1 ))
                
                if kt_chunk > length(kt_histogr[which_kt])
                    push!(kt_histogr[which_kt], zeros(Int, p.kt_nbins))
                end
                ion_velocity  = abs(x[i] - p.image_x_laser) / tof[i] / cos(3.14159 * p.measurement_angle / 180) / p.pixel_per_mm * 1e-3 # m/s
                t_surf_laser  = p.d_surf_laser * tof[i] / abs(x[i] - p.image_x_laser) * p.pixel_per_mm   # s
                offsetted_reaction_time = bl_delay[i] - t_surf_laser                                     # s
                binnr         = Int(floor(offsetted_reaction_time / p.kt_bin) + 1)
                
                if binnr > 0 && binnr <= length(kt_histogr[which_kt][kt_chunk])
                    kt_histogr[which_kt][kt_chunk][binnr] += 1 * Int64(round(ion_velocity, digits=0))
                end
            end
        end
    end
end