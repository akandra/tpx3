# 1. CREATE FOLDERS

# check if Results Folder exists, if not, generate it
if "Results" ∉ readdir(parameters.data_path) 
    mkdir(parameters.data_path * "//Results")
end

for nTOFgate in 1:size(parameters.tof_gates)[1]/2
    
    nTOFgate = Int(nTOFgate)
    
    # 2. CREATE HEADER
    header = []
    push!(header, "# --------------------------------------------------------------------------------------------\n")
    push!(header, "# experimental parameters:\n")
    for i in 1:size(exp_header)[1]
        push!(header, exp_header[i][1]*string(exp_header[i][2])*"\n")
    end
    push!(header, "# --------------------------------------------------------------------------------------------\n")
    push!(header, "# analysis parameters:\n")
    push!(header, "# TOF Gate [µs]: "*string(round.(parameters.tof_gates[nTOFgate*2-1:nTOFgate*2]*1e6, digits=2))*"\n")
    push!(header, "# t_noz_laser [µs]: "*string(round(parameters.t_noz_laser*1e6, digits=1))*"\n")
    push!(header, "# t_laser_surf [µs]: "*string(round(parameters.t_laser_surf*1e6, digits=1))*"\n")
    push!(header, "# image_x_laser : "*string(parameters.image_x_laser )*"\n")
    push!(header, "# image_y_laser : "*string(parameters.image_y_laser )*"\n")
    push!(header, "# image_x_offset: "*string(parameters.image_x_offset)*"\n")
    push!(header, "# image_y_offset: "*string(parameters.image_y_offset)*"\n")
    push!(header, "# image_x_width : "*string(parameters.image_x_width )*"\n")
    push!(header, "# image_y_width : "*string(parameters.image_y_width )*"\n")
    push!(header, "# --------------------------------------------------------------------------------------------\n")

    header_traces = copy(header)
    pop!(header_traces)
    push!(header_traces, "# mode3_t0 [s]: ",string(parameters.mode3_t0)*"\n")
    push!(header_traces, "# kt_length [s]: ",string(parameters.kt_length)*"\n")
    push!(header_traces, "# --------------------------------------------------------------------------------------------\n")
    
    # 3. SAVE KINETIC TRACE
    path2trace = parameters.data_path * "//Results//"
    open(path2trace * "Trace_TOFgate"*string(nTOFgate)*"_"*filename*".dat"; write=true) do f
        for h in header_traces
            write(f, h)
        end
        write(f, "# reaction time [ms] \t flux\n")
        writedlm(f, [x_kt reduce(hcat, traces[nTOFgate])])
    end

    # 4. SAVE TOF SPECTRUM
    open(path2trace * "TOFspectrum_" * filename * ".dat"; write=true) do f
        for h in header
            write(f, h)
        end
        write(f, "# TOS [µs] \t counts/bin\n")
        writedlm(f, [x_tof tof_spec])
    end
    
    if parameters.realtime
        # 5. SAVE REAL TIME SIGNAL
        open(path2trace * "Realtime_TOFgate"*string(nTOFgate)*"_"*filename*".dat"; write=true) do f
            for h in header
                write(f, h)
            end
            x_rt = range(parameters.realtime_t0+parameters.realtime_bin, parameters.realtime_tmax, step=parameters.realtime_bin)
            write(f, "# time [s] \t RTsig \t RTbg \n")
            writedlm(f, [x_rt rt_sig[nTOFgate] rt_bg[nTOFgate]])
        end
    end
end