using Printf
using DelimitedFiles
using Dates
using Plots
plotlyjs()

include("TPX3_analysis_structs.jl")
include("TPX3_packet_decode_functions.jl")
include("chunk_functions.jl")
include("function_rollover_correct!.jl")
include("function_analyze_kt!.jl") 
include("function_get_ROIdata.jl")
include("function_analyze_TOF!.jl")
include("function_analyze_image!.jl")
include("function_BG_histogram!.jl") ## FN ADDITIONS

function process_chunk!(BGhist_array, tof_hist_sig, tof_hist_bkg, image, kt_hist, kt_hist_bkg,
                        p, ROI_sig, ROI_bkg, 
                        out1, out2, out3, out4, out5, out6, n_out,
                        outfile_txt, ofb1, ofb2, ofb3, ofb4, ofb5, ofb6)
    
    
    println("\t beam times in this chunk: ", out1[1], " - ", out1[n_out], " seconds")
    
    if p.mode == 3 ## FN ADDITIONS
        BG_histogram!(BGhist_array, p, out1, out2, out3, out4, out5, out6, n_out, ROI_sig, ROI_bkg)
    end

    if p.create_tof
        analyze_TOF!(tof_hist_sig, p, out1, out2, out3, out4, out5, out6, n_out, ROI_sig)
        analyze_TOF!(tof_hist_bkg, p, out1, out2, out3, out4, out5, out6, n_out, ROI_bkg)
    end 

    if p.create_image
        analyze_image!(image, p, out1, out2, out3, out4, out5, out6, n_out)
    end 
    
    if p.create_kt
        analyze_kt!(kt_hist,     p, out1, out2, out3, out4, out5, out6, n_out, ROI_sig)
        analyze_kt!(kt_hist_bkg, p, out1, out2, out3, out4, out5, out6, n_out, ROI_bkg)
    end 

    write_chunk(n_out, outfile_txt, 
                        ofb1, ofb2, ofb3, ofb4, ofb5, ofb6,
                        out1, out2, out3, out4, out5, out6; 
                        write_txt=p.write_txt, write_bin = p.write_bin)
end


function convert_and_process(p :: pars)  
    
    timer = now()
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

    BGhist_array = zeros(Int32, 256) ## FN ADDITIONS
    
    # set TOFmax, TOFbin, and allocate tof histogram array
    TOFbins = Int(floor(p.tof_max/p.tof_bin))
    tof_hist_sig = zeros(Int32, TOFbins)
    tof_hist_bkg = zeros(Int32, TOFbins)

    kt_hist     = Vector{Int32}[]
    kt_hist_bkg = Vector{Int32}[]
    push!(kt_hist,     zeros(Int32, p.kt_nbins))
    push!(kt_hist_bkg, zeros(Int32, p.kt_nbins))

    image   = zeros(Int32, 256, 256)
    kt_image= zeros(Int32, 256, 256)
    kt_tof_hist= zeros(Int32, TOFbins)

    #--------------------------------------------------------------------------
    # Process data
    #--------------------------------------------------------------------------
    n_out::Int    = 0              # counter for output records
    tdc1::Float64 = 0 
    tdc2::Float64 = 0
    
# -----------------------------------------------------------------------------
#   Format of time TPX3 file:
#       Header (8 bytes) specifying size of record in bytes, chip number, mode)
#           Packet 1 (8 bytes)
#           Packet 2 (8 bytes)
#              ...
#           Packet n (8 bytes)
#       Header (8 bytes)
#           ...
#
#  Packets are of three types
#    0x4 Global time stamps, rollover time = ~ 81 days 
#    0x6 TDC data packet ,   rollover time = 107.3741824 s
#    0xB Pixel hit packet,   rollover time =  26.8435456 s
#    packet type is encoded in bits 60..63 of the packet (0 based bit indexing)
# -----------------------------------------------------------------------------

    println("------------------------------------------------------------------------------------------")
    println("Start processing file")
    println("------------------------------------------------------------------------------------------")

    # print data length in bytes
    println("filename\t\t\t: ", input_filename)
    println("packets, bytes in input file\t: ", data_length, ", ", data_length *  8)

    i = 0
    max_packet = data_length -1 # discard last packet -- it is unknown type 7 probably EOF
    
    println("testing with max_packet\t\t: ", max_packet)
    
    # loop over all the data processing one 64 bit packet word at a time 
    #       input data is handled in chunks by the getword function
    #       output data is written when an input chunk is exhausted
    get_chunk(infile, input_data)
    while i < max_packet
        i = i + 1

         ###--------------------------------------------------------------------------
            ##
            #  TODO  change to for loop
            ## 
            ###--------------------------------------------------------------------------
        
        # save and process data if we have finished processing a chunk
        if input_data.next > input_data.chunk
            
            process_chunk!(BGhist_array, tof_hist_sig, tof_hist_bkg, image, kt_hist, kt_hist_bkg,
                           p, ROI_sig, ROI_bkg, out1, out2, out3, out4, out5, out6, n_out, 
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
                # if TOA != TOA0
                #     println("correction made", TOA - TOA0)
                # end

                # save for output -- individual times
                # n_out += 1             
                # out1[n_out] = tdc1   # laser trigger
                # out2[n_out] = tdc2   # beam trigger
                # out3[n_out] = TOA    # pixel time tdc1   # laser - detection = TOF
                # out4[n_out] = TOT    # time over threshold (ns)
                # out5[n_out] = x      # pixel x
                # out6[n_out] = y      # pixel y
                
                # save for output -- differences                   
                n_out += 1             
                out1[n_out] = tdc2          # beam time
                out2[n_out] = tdc1 - tdc2   # beam to laser delay
                out3[n_out] = TOA  - tdc1   # detection time - laser = TOF
                out4[n_out] = TOT           # time over threshold (ns)
                out5[n_out] = x             # pixel x
                out6[n_out] = y             # pixel y
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
            println()
            println("---------------------------------------------")
            println("ERROR IN TPX# FILE")
            println("i = ", commas(i))
            print("packet ="); display(packet)
            println("unknown packet_type ", packet_type)
            println("---------------------------------------------")
            
            # to ignore this error, comment following line
            # stop_unknown_packet_type
        end # packet_type cascaded if-elseif selction

    end #loop over packets

    # save and process data left in chunk after reaching end of data

    if n_out > 0
        process_chunk!(BGhist_array, tof_hist_sig, tof_hist_bkg, image, kt_hist, kt_hist_bkg,
                       p, ROI_sig, ROI_bkg, out1, out2, out3, out4, out5, out6, n_out, 
                       outfile_txt, ofb1, ofb2, ofb3, ofb4, ofb5, ofb6)
    end
                
    println()
    println("done converting data\n")
    println("------------------------------------------------------------------------------------------")
    println("Processing time\t\t\t: ", now() - timer)
    println("Total exp. time\t\t\t: ", out1[n_out], " seconds")
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

    if p.mode == 3
        println("\n Plotting Background Histogram for the first seconds")  ## FN ADDITIONS
        vlines = [p.image_x_laser - p.image_x_offset - p.image_x_width,
                    p.image_x_laser - p.image_x_offset,
                    p.image_x_laser,
                    p.image_x_laser + p.image_x_offset,
                    p.image_x_laser + p.image_x_offset + p.image_x_width]
        
        p0 = plot(range(1, 256, step=1), BGhist_array)
        p0 = vline!(transpose(vlines), linecolor=[:red :red :black :green :green])
        display(p0)
    end               
                
                  
    # results of TOF analysis
    # println("\nrange= ", range(0, TOFmax, length(tof_hist)))
    # println("TOFmax = ", TOFmax, "  length = ", length(tof_hist))
    TOF = range(0, p.tof_max * 1e6, TOFbins)

    # create box to show TOF gate
    imin = Int(floor(p.tof_gate_min/p.tof_bin))
    imax = Int(floor(p.tof_gate_max/p.tof_bin))
    TOFgate_max = maximum(tof_hist_sig[imin:imax])
    TOFgate_box = [p.tof_gate_min p.tof_gate_min  p.tof_gate_max p.tof_gate_max
                   0 TOFgate_max TOFgate_max 0]

    # plot complete mass spectrum
    println("Plotting mass spectrum")
    p1 = plot( TOF, [-tof_hist_sig -tof_hist_bkg (tof_hist_sig - tof_hist_bkg)], 
                xlabel="TOF (µs)",ylabel="counts / bin",
                labels = ["-sig" "-bkg" "sig-bkg"], framestyle = :box,
                title = "TOF mass spectrum")
    p1 = plot!(p1, TOFgate_box[1,:]*1e6, TOFgate_box[2,:], labels= :none, ann=( (0.99,0.05), text(p.filename_stem,:right,8)))
    display(p1)

    # # plot gates part of mass spectrum (+/-200ns)
    # plot( TOF, tof_hist, 
    #             xlabel="TOF (µs)",ylabel="counts / bin",
    #             xlim = (TOFgate[1]*1e6 - 0.2, TOFgate[2]*1e6 + 0.2),
    #             ylim = (0, TOFgate_max *1.2),
    #             legend = :none, framestyle = :box,
    #             title = "TOF mass spectrum" )          
    # display(plot!(TOFgate_box[1,:]*1e6, TOFgate_box[2,:]))

    # plot image
    println("\n Sum Image for TOF Gate")
    p1 = Plots.heatmap(image, clims=p.image_scale, ann=( (0.95,0.95), text(p.filename_stem,:white,:right,8)))
    
    plot!(p1,  [ROI_sig[1,1], ROI_sig[1,2], ROI_sig[1,2], ROI_sig[1,1], ROI_sig[1,1]],
               [ROI_sig[2,1], ROI_sig[2,1], ROI_sig[2,2], ROI_sig[2,2], ROI_sig[2,1]],
               color=:green)
    plot!(p1,  [ROI_bkg[1,1], ROI_bkg[1,2], ROI_bkg[1,2], ROI_bkg[1,1], ROI_bkg[1,1]],
               [ROI_bkg[2,1], ROI_bkg[2,1], ROI_bkg[2,2], ROI_bkg[2,2], ROI_bkg[2,1]],
                color=:red)
    display(p1)

    #--------------------------------------------------------------------------
    # Plot the Kinetic trace
    #--------------------------------------------------------------------------
    if p.create_kt    
        println("\nPlotting Kinetic Traces in Chunks of ", p.kt_length, " seconds")
        # subtracting background from the kinetic trace
        kt_hist_BGcorr = kt_hist .- kt_hist_bkg

        # define the time range
        t = range(0, p.kt_max * 1e3, p.kt_nbins)

        kt_norm = push!([ p.kt_length for i in 1:length(kt_hist)-1], mod(out1[n_out] - p.kt_t0, p.kt_length))
        #println("\t\tNormalization vector ", kt_norm)
        
        endvals   = cumsum(kt_norm) .+ parameters.kt_t0
        startvals = endvals .- kt_norm

        display(plot(t , kt_hist_BGcorr ./ kt_norm,
            xlabel="reaction time (ms)",ylabel="counts / bin / s",
#            labels = reshape(legend, 1, size(legend)[1]),
            framestyle = :box,
            title = "Kinetic traces",
            ann=( (0.99,0.95), text(p.filename_stem,:right,8))
            ))

    end
end # function read_and_convert