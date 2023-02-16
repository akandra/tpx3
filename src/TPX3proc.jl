using Parameters

include("function_convert_and_process.jl")
include("TPX3_analysis_structs.jl")
version = "0.1.9"

# define parameters

parameters = pars()

parameters.data_path      = "../Data/20221010" # select the data file to process
data_filenames            = readdir(parameters.data_path)
TPX3_fnames               = filter(x-> occursin("tpx3",x) ,data_filenames)

fn_n = 6

filename                  = split(TPX3_fnames[fn_n], '.')[1]
parameters.filename_stem  = parameters.data_path * "/" * filename

parameters.write_bin      = false
parameters.write_txt      = false

parameters.mode           = 3 ## FN ADDITIONS, has to be completed

# TOF parameters
parameters.create_tof = true
parameters.tof_max = 8000e-9
parameters.tof_bin = 20e-9
parameters.tof_gate_min = 1.0e-6#5630e-9
parameters.tof_gate_max = 8000.0e-9#5900e-9

# image parameters
parameters.create_image  = true
parameters.image_scale = (0,20)
parameters.image_x_laser  = 140
parameters.image_y_laser  = 172
parameters.image_x_offset = 5
parameters.image_y_offset = 0
parameters.image_x_width = 25
parameters.image_y_width = 20

# Kinetic Trace create parameters

parameters.create_kt     = true
parameters.kt_max = 2e-3
parameters.kt_bin = 20e-6
parameters.kt_nbins = Int(floor(parameters.kt_max/parameters.kt_bin))
parameters.kt_gate_min = 0e-6
parameters.kt_gate_max = 100000e-6
parameters.kt_t0 = 5.0
parameters.kt_length = 110.0

# Laser and nozzle frequencies (Hz)
parameters.freq_nozzle = 20.
parameters.freq_laser  = 100000.

parameters.chunk_size     = (500*1024*1024)÷8

println("selected file\t\t: ", filename)

convert_and_process(parameters)

#NEXTTIME
# 1. time-dependent kinetics - finish the normalization and the legend
# 1. deal with border case fo get_chunk when there is no leftover
# 2. Multimass KT
# 3. Introduce the workflow keys: 
#        1: calibration of tof from background;
#        2: getting incident beam profile and select mass
#        3: getting kinetic traces

# Long-term perspective:
# 1. julia package
