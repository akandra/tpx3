# tpx3 file naming convention:
#[file number]_[metal]_[facet]_[beam species]_[beam RR]_[laser RR]_[laser position]_[temperature]_[comment]_xxxxxx.tpx3

#using Gtk
using Parameters

include("function_convert_and_process.jl")
include("TPX3_analysis_structs.jl")
version = "0.1.10"

# define parameters

parameters = pars()

# PATH 1
parameters.data_path      = "F://B2Data//09.02.2023"
filename                  = "004_RealTimeTitrations_Ptlower_480C_50HzNH3_saturationO2_10minDosingAt100Hz_000000"

# PATH 2
parameters.data_path      = "/home/akandra/Dropbox/Timepix camera data analysis/Data"
filename                  = "20221010/004_100kHz_Ptlowerfacet_thickNOat20Hz_480C000000"

parameters.filename_stem  = parameters.data_path * "/" * filename

#filename = nothing
#if isnothing(filename)
#    filename = open_dialog("Select a TPX3 file")
#    parameters.filename_stem  = split(filename, ".TPX3")[1]
#end


#-------------------------------------------------------------------#
#-----------------        PARAMETER SECTION        -----------------#

parameters.write_bin      = false
parameters.write_txt      = false

parameters.mode           = 3
parameters.first_seconds  = 3.0

# TOF parameters
parameters.tof_max        = 8e-6
parameters.tof_bin        = 10e-9
parameters.tof_gates      = [ 5.63, 5.90, 5.95, 6.0 ] .* 1e-6


# image parameters 
parameters.create_image   = true

parameters.image_x_laser  = 140
parameters.image_y_laser  = 172
parameters.image_x_offset = 5
parameters.image_y_offset = 0
parameters.image_x_width  = 25
parameters.image_y_width  = 20

# Kinetic Trace create parameters
parameters.kt_max         = 2e-3
parameters.kt_bin         = 20e-6
parameters.kt_nbins       = Int(floor(parameters.kt_max/parameters.kt_bin))
parameters.kt_t0          = 5.0
parameters.kt_length      = 10.0

# Laser and nozzle frequencies (Hz)
parameters.freq_nozzle    = 20.
parameters.freq_laser     = 100000.

parameters.chunk_size     = (500*1024*1024)÷8

println("selected file: ", filename)

convert_and_process(parameters)

#NEXTTIME
# 1. use a profiler to improve the performance
# 2. deal with border case fo get_chunk when there is no leftover
# 3. Introduce more modes
#       a) imitate on-the-fly processing of the data
# Long-term perspective:
# 1. julia package
