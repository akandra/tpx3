using Printf
using DelimitedFiles
using Dates
using Plots
plotlyjs()

include("TPX3_analysis_structs.jl")
include("TPX3_packet_decode_functions.jl")
include("chunk_functions.jl")
include("function_rollover_correct!.jl")
include("function_get_ROIdata.jl")
include("function_analysis_mode2!.jl")
include("function_analysis_mode3!.jl")
include("function_BG_histogram!.jl")
include("function_analyze_image!.jl")
include("function_realtime_analysis!.jl")

function process_chunk!(BGhists, tof_hist_sig, tof_hist_bkg, images, 
                        kt_hist_sig,kt_hist_bkg,
                        realtime_trace_sig, realtime_trace_bkg,
                        kt_firstsec_sig, kt_firstsec_bkg, p, ROI_sig, ROI_bkg, ROI_beam, 
                        out1, out2, out3, out4, out5, out6, n_out,
                        outfile_txt, ofb1, ofb2, ofb3, ofb4, ofb5, ofb6)
    
    
    println("\t beam times in this chunk: ", out1[1], " - ", out1[n_out], " seconds")
    
    
    # MODE 3 (KINETIC TRACE) ANALYSIS
    
    if p.mode == 3
        if out1[1] < p.mode3_tmax
            BG_histogram!(BGhists, p, out1, out2, out3, out4, out5, out6, n_out, ROI_sig, ROI_bkg)
            analysis_mode3!(tof_hist_sig, kt_hist_sig, kt_firstsec_sig, p, ROI_sig, out1, out2, out3, out4, out5, out6, n_out)
            analysis_mode3!(tof_hist_bkg, kt_hist_bkg, kt_firstsec_bkg, p, ROI_bkg, out1, out2, out3, out4, out5, out6, n_out)

            if p.create_image
                analyze_image!(images, p, out1, out2, out3, out4, out5, out6, n_out)
            end
        end
    end
    
    # REAL TIME ANALYSIS
    
    if p.realtime
        if out1[1] < p.realtime_tmax || out1[n_out] > p.realtime_t0
            realtime_analysis!(realtime_trace_sig, p, ROI_sig, out1, out2, out3, out4, out5, out6, n_out)
            realtime_analysis!(realtime_trace_bkg, p, ROI_bkg, out1, out2, out3, out4, out5, out6, n_out)
        end
    end
    
    
    # MODE 2 (INCIDENT BEAM) ANALYSIS
    
    if p.mode == 2
        analysis_mode2!(tof_hist_sig, kt_hist_sig, kt_firstsec_sig, p, ROI_beam, out1, out2, out3, out4, out5, out6, n_out)

        if p.create_image
            analyze_image!(images, p, out1, out2, out3, out4, out5, out6, n_out)
        end
        
    end
    
    # WRITING OUTPUT FILES
    
    write_chunk(n_out, outfile_txt, 
                        ofb1, ofb2, ofb3, ofb4, ofb5, ofb6,
                        out1, out2, out3, out4, out5, out6; 
                        write_txt=p.write_txt, write_bin = p.write_bin)
end


function convert_and_process(p :: pars)  
    
    timer = time()
    # signature of record header for checking
    hdr_check::UInt32 = reinterpret(UInt32, Vector{UInt8}("TPX3"))[1]
    

    input_filename        = p.filename_stem * ".tpx3"
    output_filename_txt   = p.filename_stem        * "_jl_" * version * ".txt"
    output_filename_bin1  = p.filename_stem * "_1" * "_jl_" * version * ".bin"
    output_filename_bin2  = p.filename_stem * "_2" * "_jl_" * version * ".bin"
    output_filename_bin3  = p.filename_stem * "_3" * "_jl_" * version * ".bin"
    output_filename_bin4  = p.filename_stem * "_4" * "_jl_" * version * ".bin"
    output_filename_bin5  = p.filename_stem * "_5" * "_jl_" * version * ".bin"
    output_filename_bin6  = p.filename_stem * "_6" * "_jl_" * version * ".bin"
    
    input_data = get_word_state(p.chunk_size, 0, 1, zeros(UInt64, p.chunk_size));

    infile = open(input_filename)
    seekend(infile)
    # get length of the data in 8 byte words
    data_length = position(infile) ÷ 8
    seekstart(infile)

    # open files and initialize timers
    outfile_txt = p.write_txt ? open(output_filename_txt, "w") : nothing
    ofb1 = p.write_bin ? open(output_filename_bin1, "w") : nothing
    ofb2 = p.write_bin ? open(output_filename_bin2, "w") : nothing
    ofb3 = p.write_bin ? open(output_filename_bin3, "w") : nothing
    ofb4 = p.write_bin ? open(output_filename_bin4, "w") : nothing
    ofb5 = p.write_bin ? open(output_filename_bin5, "w") : nothing
    ofb6 = p.write_bin ? open(output_filename_bin6, "w") : nothing

    # initialize rollover correction related variable
    # values of previous time
    TDC_prior  = 0.0
    TOA_prior  = 0.0
    
    # corrections to add
    TDC_add  = 0.0
    TOA_add  = 0.0

    # rollover amounts
    TDC_tick   = 1/320e6  # time bin from 320 MHz.  clock
    TDC_ro     = 2^35 * TDC_tick
    TOA_tick   = 25e-9
    TOA_ro     = 2^30 * TOA_tick
    
    #--------------------------------------------------------------------------
    # allocate vectors to store output results  
    out1 = zeros(Float32, p.chunk_size)
    out2 = zeros(Float32, p.chunk_size)
    out3 = zeros(Float32, p.chunk_size)
    out4 = zeros(UInt16 , p.chunk_size)
    out5 = zeros(UInt8  , p.chunk_size)
    out6 = zeros(UInt8  , p.chunk_size)

    # define ROIs for the signal and background

    laser_pos = (p.image_x_laser,  p.image_y_laser)
    offset    = (p.image_x_offset, p.image_y_offset)
    width     = (p.image_x_width,  p.image_y_width)

    ROI_sig = [ [laser_pos[1] + offset[1] laser_pos[1] + offset[1] + width[1] ]; 
                [laser_pos[2] + offset[2] - width[2] laser_pos[2] + offset[2] + width[2] ]]

    ROI_bkg = [ [laser_pos[1] - offset[1] - width[1] laser_pos[1] - offset[1] ]; 
                [laser_pos[2] + offset[2] - width[2] laser_pos[2] + offset[2] + width[2] ]]

    # ROI for incident beam
    
    ROI_beam  = [[p.beam_x_min p.beam_x_max]; [p.beam_y_min p.beam_y_max]]
    
    ##################   INITIALIZE ARRAYS   ##################
    
    TOFbins            = Int(floor(p.tof_max/p.tof_bin))
    tof_hist_sig       = zeros(Int32, TOFbins)
    tof_hist_bkg       = zeros(Int32, TOFbins)

    nr_realtimebins    = Int(floor((p.realtime_tmax - p.realtime_t0)/p.realtime_bin))        
    
    n_tof_gates = Int(size(p.tof_gates)[1]/2)
    # if number of kinetic traces  > 20, then push! is used to increase the size of array
    kt_hist_sig        = [[ zeros(Int, p.kt_nbins) for _ in 1:20] for _ in 1:n_tof_gates]
    kt_hist_bkg        = [[ zeros(Int, p.kt_nbins) for _ in 1:20] for _ in 1:n_tof_gates]

    kt_firstsec_sig    = [ zeros(Int,  p.kt_nbins) for _ in 1:n_tof_gates]
    kt_firstsec_bkg    = [ zeros(Int,  p.kt_nbins) for _ in 1:n_tof_gates]
    
    realtime_trace_sig = [ zeros(Int,  nr_realtimebins) for _ in 1:n_tof_gates]
    realtime_trace_bkg = [ zeros(Int,  nr_realtimebins) for _ in 1:n_tof_gates]

    BGhists            = [ zeros(Int, 256)      for _ in 1:n_tof_gates]
    images             = [ zeros(Int, 256, 256) for _ in 1:n_tof_gates]
 
    ############################################################
    
    #--------------------------------------------------------------------------
    # Process data
    #--------------------------------------------------------------------------
    n_out::Int    = 0              # counter for output records
    tdc1::Float64 = 0 
    tdc2::Float64 = 0
    println("filesize\t: ", round(filesize(input_filename) * 9.313e-10, digits=3), " GB")
    println("------------------------------------------------------------------------------------------")
    println("Start processing file")
    println("------------------------------------------------------------------------------------------")

    i = 0
    max_packet = data_length -1 # discard last packet -- it is unknown type 7 probably EOF
    
    
    # Introduce a variable to stop the loading/converting procedure if one is not interested in the whole file
    early_stop = 0
    
    get_chunk(infile, input_data)
    
    while i < max_packet && early_stop < 1
        i = i + 1
        
        
        # -------------------------------------------------- #
        # COMMENTED OUT FOR COMBINED RRR+TITRATION EXPERIMENTS
        
        #if out1[1] > p.mode3_tmax
        #    early_stop += 1
        #end
        # -------------------------------------------------- #
        
        # save and process data if we have finished processing a chunk
        if input_data.next > input_data.chunk
            
           process_chunk!(BGhists, tof_hist_sig, tof_hist_bkg, images, 
                        kt_hist_sig, kt_hist_bkg, realtime_trace_sig, realtime_trace_bkg,
                        kt_firstsec_sig, kt_firstsec_bkg, p, ROI_sig, ROI_bkg, ROI_beam, 
                        out1, out2, out3, out4, out5, out6, n_out,
                        outfile_txt, ofb1, ofb2, ofb3, ofb4, ofb5, ofb6)
            n_out = 0

            # get a new chunk
            get_chunk(infile, input_data)
            ###--------------------------------------------------------------------------
            ##
            #   TODO  figure out if we need an empty chunk test
            #         Probably, we need
            ## 
            ###--------------------------------------------------------------------------
        end

        packet = input_data.buffer[input_data.next]
        input_data.next  += 1
        
        # ------------------------------------------------------------------------
        #   get the packet type
        #   packet type is in bits 60 to 63 of data packet (bits numbered 0 .. 63)
        # ------------------------------------------------------------------------
        packet_type = packet >> 60
        
        #-------------------------------------------------------------------------
        #   Check the hdear using signature of header
        # ------------------------------------------------------------------------
        if UInt32(packet & 0xffffffff) == hdr_check
            chipnr = (packet >> (4*8)) & 0xff
            mode   = (packet >> (5*8)) & 0xff
            size   = (packet >> (6*8)) & 0xffff
            
            # number of packets for this header is size(bytes) / 8
            n_pack = size >> 3 
            
        #------------------------------------------------------------------
        # check for trigger timestamp packet: extract tdc1 and tdc2
        #------------------------------------------------------------------
        elseif packet_type == 0x6
            type, time = process_TDC_timestamp_packet(packet)
            time, TDC_prior, TDC_add = rollover_correct!(time, TDC_prior, TDC_ro, TDC_add)
                
            # save timestamp for output
            if type == 0x6f
                tdc1 = time
            elseif type == 0x6e 
                tdc2 = time
            end
            
        #------------------------------------------------------------------
        # check for chip data packet: extract ToA and ToT timestamp, x, y
        #------------------------------------------------------------------
        elseif packet_type == 0xb
            # ignore data if no molecular beam trigger received
            if tdc2 > 0.0
                TOA, TOT, x, y = process_pixel_hit_packet(packet)
                TOA0 = TOA
                TOA, TOA_prior, TOA_add = rollover_correct!(TOA, TOA_prior, TOA_ro, TOA_add)

                # save for output -- differences                   
                n_out += 1             
                out1[n_out] = tdc2              # beam time
                out2[n_out] = tdc1 - tdc2       # beam to laser delay
                out3[n_out] = TOA  - tdc1       # detection time - laser = TOF
                out4[n_out] = TOT               # time over threshold (ns)
                out5[n_out] = x                 # pixel x
                out6[n_out] = y                 # pixel y
            end
            
        #------------------------------------------------------------------
        # check for global timestamp packet
        #------------------------------------------------------------------
        elseif packet_type == 0x4
            
            process_global_timestamp_packet(packet)
            # when we figure out the global timestamp we have add it to the ouput
        #------------------------------------------------------------------
        # unknown packet type --> error
        # ASI c++ code ignores this error
        #------------------------------------------------------------------
        else
            #println()
            #println("---------------------------------------------")
            #println("ERROR IN TPX# FILE")
            #println("i = ", commas(i))
            #print("packet ="); display(packet)
            #println("unknown packet_type ", packet_type)
            #println("---------------------------------------------")
            
            # to ignore this error, comment following line
            # stop_unknown_packet_type
        end # packet_type cascaded if-elseif selction

    end #loop over packets
    # save and process data left in chunk after reaching end of data

    if n_out > 0
       process_chunk!(BGhists, tof_hist_sig, tof_hist_bkg, images, 
                      kt_hist_sig, kt_hist_bkg, realtime_trace_sig, realtime_trace_bkg,
                      kt_firstsec_sig, kt_firstsec_bkg, p, ROI_sig, ROI_bkg, ROI_beam,
                      out1, out2, out3, out4, out5, out6, n_out,
                      outfile_txt, ofb1, ofb2, ofb3, ofb4, ofb5, ofb6)
    end
                
    println()
    println("done converting data\n")
    println("------------------------------------------------------------------------------------------")
    total_proc_time = time() - timer
    filesize_GB     = filesize(input_filename) * 9.313e-10
    println("Processing time")
    println("\tAbsolute: ", round(total_proc_time, digits=1), "\tseconds")
    println("\tRelative: ", round(total_proc_time / 60 / filesize_GB, digits=2), "\tmin/GB")
    println("Total exp. time\t: ", round(out1[n_out], digits=1), "\tseconds")
    println("------------------------------------------------------------------------------------------")
                
    #--------------------------------------------------------------------------
    #   Done converting file
    #--------------------------------------------------------------------------   
    close(infile)
    if p.write_txt
        close(outfile_txt)
    end
    if p.write_bin
        close(ofb1)
        close(ofb2)
        close(ofb3)
        close(ofb4)
        close(ofb5)
        close(ofb6)
    end

                
    # --------------------------------------------------------------------- #
    # -------------------------  PLOTTING SECTION ------------------------- #
                
                
    ###########   PLOTTING BACKGROUND FOR FIRST SECONDS   ##########
                
    if p.mode == 3
        legend_nr = 1
        for nr in 1:Int(size(p.tof_gates)[1]/2)
            vlines = [p.image_x_laser - p.image_x_offset - p.image_x_width,
                        p.image_x_laser - p.image_x_offset,
                        p.image_x_laser,
                        p.image_x_laser + p.image_x_offset,
                        p.image_x_laser + p.image_x_offset + p.image_x_width]

            p0 = plot(range(1, 256, step=1), BGhists[nr], label="TOF Gate " * string(round(p.tof_gates[legend_nr]*1e6, digits=2)) 
                                * " - " * string(round(p.tof_gates[legend_nr+1]*1e6, digits=2)) * " μs")
            p0 = vline!(transpose(vlines), linecolor=[:red :red :black :green :green], labels=:none)
            p0 = title!("Background Histogram of first " * string(round(p.first_seconds, digits=1)) * " seconds")
            p0 = plot!(legend=:right)
            display(p0)
            legend_nr += 2
        end  
    end

                
    ###################   PLOTTING TOF SPECTRUM   ##################
                
    if p.mode == 2 || p.mode == 3         
        xtof    = range(0, p.tof_max * 1e6, TOFbins)
        tofspec = tof_hist_sig .- tof_hist_bkg

        p1      = plot(xtof, tofspec, label=:none)

        for nr in 1:size(p.tof_gates)[1]
            if mod(nr, 2) == 1
                tofbox = [p.tof_gates[nr] p.tof_gates[nr] p.tof_gates[nr+1] p.tof_gates[nr+1] 
                0 findmax(tofspec)[1] findmax(tofspec)[1] 0]
                p1 = plot!(tofbox[1,:]*1e6, tofbox[2,:], labels=:none)
            end
        end
        if p.mode == 3
            p1 = title!("background corrected TOF Spectrum")
        end
        if p.mode == 2
            p1 = title!("TOF Spectrum")
        end
        p1 = xlabel!("TOF / μs") 
        p1 = ylabel!("counts / bin")         
        display(p1)     
    end            
                
    
    ######################   PLOTTING SUM IMAGES   #####################
                
    if p.create_image  
        legend_nr = 1 
        for nr in 1:Int(size(p.tof_gates)[1]/2)
            p2 = plot()
            p2 = title!("Sum Image in TOF Gate: " * string(round(p.tof_gates[legend_nr]*1e6, digits=2)) * " - " 
                    * string(round(p.tof_gates[legend_nr+1]*1e6, digits=2)) * " μs")
            
            if p.mode == 3
            p2 = heatmap!(images[Int(nr)], labels=:none)
            plot!(p2,  [ROI_sig[1,1], ROI_sig[1,2], ROI_sig[1,2], ROI_sig[1,1], ROI_sig[1,1]],
                   [ROI_sig[2,1], ROI_sig[2,1], ROI_sig[2,2], ROI_sig[2,2], ROI_sig[2,1]],
                   color=:green, labels=:none)
            plot!(p2,  [ROI_bkg[1,1], ROI_bkg[1,2], ROI_bkg[1,2], ROI_bkg[1,1], ROI_bkg[1,1]],
                       [ROI_bkg[2,1], ROI_bkg[2,1], ROI_bkg[2,2], ROI_bkg[2,2], ROI_bkg[2,1]],
                        color=:red, labels=:none)
            plot!(p2, [p.image_x_laser, p.image_x_laser], [0, 256], color=:yellow, labels=:none)
            display(p2)
            end
                        
            
            if p.mode == 2
            plot!(p2,  [ROI_beam[1,1], ROI_beam[1,2], ROI_beam[1,2], ROI_beam[1,1], ROI_beam[1,1]],
                       [ROI_beam[2,1], ROI_beam[2,1], ROI_beam[2,2], ROI_beam[2,2], ROI_beam[2,1]],
                        color=:red, labels=:none)
            plot!(p2, [p.image_x_laser, p.image_x_laser], [0, 256], color=:yellow, labels=:none)
            p2 = heatmap!(images[Int(nr)], labels=:none, clim=(0, p.beam_zscale))
            display(p2)
            legend_nr += 2
                        
            DeltaX    = abs(ROI_beam[1,1] - ROI_beam[1,2])
            DeltaTOF  = abs(p.tof_gates[nr+1] - p.tof_gates[nr])
            TOFmean   = p.tof_gates[nr] + DeltaTOF/2
            Xmean     = p.image_x_laser - ROI_beam[1,2] + DeltaX/2
            vbeam     = Xmean/cos(3.141592 * p.sidenoz_angle / 180)/TOFmean/p.pixel_per_mm*1e-3
            rel_error = DeltaX/2/Xmean + DeltaTOF/2/TOFmean
            println("incident beam velocity         = (", Int(round(vbeam, digits=0)), 
                            " +- ", Int(round(vbeam*rel_error, digits=0)), ") m/s")
                        
            flight_time_us = p.d_surf_laser/cos(3.141592 * p.sidenoz_angle / 180) * 1e-3 / vbeam * 1e6
                        
            println("flight time (laser to surface) = (", Int(round(flight_time_us, digits=0)), " +- ", 
                            Int(round(flight_time_us*rel_error, digits=0)), ") us")
            end
        end
    end
            
   # CORRECTION: Normalization Introduced on 12.04.2024
   kt_hist_sig        = kt_hist_sig ./ p.kt_bin
   kt_hist_bkg        = kt_hist_bkg ./ p.kt_bin        
   realtime_trace_sig = realtime_trace_sig ./ p.realtime_bin
   realtime_trace_bkg = realtime_trace_bkg ./ p.realtime_bin
                            
    ###################   PLOTTING KINETIC TRACES   ################### 
                
    if p.mode == 2 || p.mode == 3 
        endtime     = findmin([p.mode3_tmax, out1[n_out]])[1]
        nr_traces   = ceil((endtime-p.mode3_t0) / p.kt_length)                   
        kt_hist_sig = [kt_hist_sig[i][1:Int(nr_traces)] for i in 1:size(kt_hist_sig)[1]]
        kt_hist_bkg = [kt_hist_bkg[i][1:Int(nr_traces)] for i in 1:size(kt_hist_bkg)[1]]
        traces      = kt_hist_sig .- kt_hist_bkg    
                       
        # Calculate a normalization vector for the kinetic traces (24.05.2024)
        #  traces should be normalized to the number of beam shots for which the flux per kinetic trace was summed up                      
        kt_maxtime  = findmin([out1[n_out], p.mode3_tmax])[1]
        kt_times    = push!([p.kt_length for i in 1:length(traces[1])-1], mod(kt_maxtime - p.mode3_t0, p.kt_length))
        norm_vec    = kt_times .* p.freq_nozzle
                    
        endvals     = cumsum(kt_times) .+ p.mode3_t0
        startvals   = endvals .- kt_times
        legend      = [string(round(startvals[i],digits=1))*" - "*string(round(endvals[i],digits=1))*" s" for i in 1:size(startvals)[1]]
             
        if p.mode == 2
            x_kt   = range(0, p.kt_max * 1e3, p.kt_nbins)
        end
        if p.mode == 3
           # Now correct for the reaction time
           x_kt    = range(0, p.kt_max * 1e3, p.kt_nbins) .- (p.t_noz_laser + p.t_laser_surf) * 1e3
        end
                    
        legend_nr = 1
        for j in 1:size(traces)[1]
            p3 = plot()
            p3 = plot!(x_kt, traces[j] ./ norm_vec, labels = reshape(legend, 1, size(legend)[1]))
            p3 = title!("TOF Gate: " * string(round(p.tof_gates[legend_nr]*1e6, digits=2)) * " -> " 
                * string(round(p.tof_gates[legend_nr+1]*1e6, digits=2)) * " μs")
                        
            if p.mode == 2
                p3 = xlabel!("beam laser delay / ms")
                p3 = ylabel!("counts")
            end  
            if p.mode == 3
                p3 = xlabel!("reaction time / ms")
                p3 = ylabel!("flux")
            end      
            display(p3)
            legend_nr += 2
        end          
    end 
   
# units of traces:         flux / s (/ beam shot)
# units of realtime_trace: flux / s
   
   return BGhists, xtof, tofspec, x_kt, [traces[i] ./ norm_vec for i in 1:size(traces)[1]], [kt_hist_bkg[i] ./ norm_vec for i in 1:size(kt_hist_bkg)[1]], realtime_trace_sig, realtime_trace_bkg, out1[n_out] # , (kt_firstsec_sig .- kt_firstsec_bkg) ./ p.first_seconds,                 
end