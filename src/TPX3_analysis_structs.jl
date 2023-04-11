# Program units:
# time is in seconds
# frequiency is in 1/seconds

mutable struct get_word_state
    chunk       :: Int              # = 1024*1024*1000
    chunk_number:: Int              # = 0
#    first_call  :: Bool             # = true
    next        :: Int              # = 0
    buffer      :: Vector{UInt64}   # = zeros(chunk)
end

@with_kw mutable struct pars

    data_path       :: String = "./"
    filename_stem   :: String = ""

    write_bin       :: Bool = false
    write_txt       :: Bool = false
    
    mode            :: Int16   = 0
    mode3_t0        :: Float64 = 0
    mode3_tmax      :: Float64 = Inf

    tof_max         :: Float64 = 0.0
    tof_bin         :: Float64 = 0.0
    
    tof_gates       :: Vector{Float64} = []
    first_seconds   :: Float64 = 3.0
    
    create_image    :: Bool  = false   
    image_x_laser   :: Int16 = 0
    image_y_laser   :: Int16 = 0
    image_x_offset  :: Int16 = 0
    image_y_offset  :: Int16 = 0
    image_x_width   :: Int16 = 0
    image_y_width   :: Int16 = 0
   
    beam_x_min     :: Int16 = 0
    beam_x_max     :: Int16 = 0
    beam_y_min     :: Int16 = 0
    beam_y_max     :: Int16 = 0
    beam_zscale    :: Int16 = 0

    kt_max          :: Float64 = 0.0
    kt_bin          :: Float64 = 0.0
    kt_nbins        :: Int = 0
    kt_gate_min     :: Float64 = 0.0
    kt_gate_max     :: Float64 = Inf
    kt_length       :: Float64 = Inf
    t_noz_laser     :: Float64 = 0
    t_laser_surf    :: Float64 = 0
    

    realtime        :: Bool    = false
    realtime_t0     :: Float64 = 0.0
    realtime_bin    :: Float64 = 0.0
    realtime_tmax   :: Float64 = Inf
    
    freq_nozzle     :: Float64 = 20.
    freq_laser      :: Float64 = 100000.
    
    # Geometries of the Beamer3 Machine
    pixel_per_mm      :: Float64 = 2 * 103/49  # pixel / mm
    sidenoz_angle     :: Float64 = 30          # angle between sidenozzle and surface normal (degree)
    d_surf_laser      :: Float64 = 25          # mm
    measurement_angle :: Float64 = 0           # angle between measurement position and surface normal (degree)

    chunk_size      :: Int = (1024*1024*2000)÷8
end
