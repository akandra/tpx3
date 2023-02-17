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
    
    mode            :: Int16 = 0

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

    kt_max          :: Float64 = 0.0
    kt_bin          :: Float64 = 0.0
    kt_nbins        :: Int = 0
    kt_gate_min     :: Float64 = 0.0
    kt_gate_max     :: Float64 = Inf
    kt_t0           :: Float64 = 0.0
    kt_length       :: Float64 = Inf

    freq_nozzle     :: Float64 = 20.
    freq_laser      :: Float64 = 100000.

    chunk_size      :: Int = (1024*1024*2000)÷8
end
