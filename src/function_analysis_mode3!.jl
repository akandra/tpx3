function find_TOF_interval!(tof_value, p)
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



function analysis_mode3!(tof_hist, kt_histogr, kt_histogr_first_sec, p, ROI, beam_time, bl_delay, tof, tot, x, y, n_out)
    
    idx_roi       =  get_ROIdata(beam_time, bl_delay, tof, tot, x, y, n_out, ROI)

    # NOW, LOOP ONCE OVER THE EVENTS IN THE ROI
    for i in idx_roi
        
        # 1. GENERATE TOF SPECTRUM
        tofidx   =   Int(floor(tof[i] / p.tof_bin) + 1)    
        if tofidx > 1  &&  tofidx <= length(tof_hist)
            tof_hist[tofidx] += 1
        end
        
        # 2. GENERATE KINETIC TRACE FOR FIRST SECONDS
        if beam_time[i] > 0 && beam_time[i] < p.first_seconds
            
            which_kt = find_TOF_interval!(tof[i], p)
            
            if typeof(which_kt) != Nothing
                binnr = Int(floor(bl_delay[i] / p.kt_bin) + 1)
                
                if binnr > 0 && binnr <= length(kt_histogr_first_sec[which_kt])
                    kt_histogr_first_sec[which_kt][binnr] += 1
                end
            end
        end
        
        # 3. GENERATE KINETIC TRACES IN CHUNKS FOR TIMES > kt_t0
        if beam_time[i] > p.kt_t0
            
            which_kt = find_TOF_interval!(tof[i], p)
            
            if typeof(which_kt) != Nothing
                kt_chunk = Int(floor( (beam_time[i] - p.kt_t0) / p.kt_length + 1 ))
                
                if kt_chunk > length(kt_histogr[which_kt])
                    push!(kt_histogr[which_kt], zeros(Int32, p.kt_nbins))
                end
                
                binnr = Int(floor(bl_delay[i] / p.kt_bin) + 1)
                
                if binnr > 0 && binnr <= length(kt_histogr[which_kt][kt_chunk])
                    kt_histogr[which_kt][kt_chunk][binnr] += 1
                end
            end
        end
    end
end          