"""
add commas to string representing an integer
"""
function commas(num::Integer)
    str = string(num)
    return replace(str, r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
end




"""
Process tdc timestamp packet 

return type, time


"""
# types: ``0x6f → 1+, 0x6a \→ 1-, 0x6e \→ 2+, 0x6b \→ 2-``


function process_TDC_timestamp_packet(packet)
    
    # get the trigger type from bits 57-64 (1 based)
    type = packet >> 56   
        
    time = Float64( (packet >> 9 & 0x7FFFFFFFF)) * 3.125E-9 +   #coarse_time
           Float64( (packet >> 5 & 0xf) - 1)   * 260E-12        #fine_time
    
    return type, time
end # function process_TDC_timestamp_packet    





"""
Process Chip data packet

return ToA, ToT, x, y
"""
function process_pixel_hit_packet(packet)

    TOA  = Float64( packet >> 30 & 0x3fff  +    # TOA field
                   (packet & 0xffff) * 16384    # spidrTime field
                  ) * 25E-9 -                    # convert to ns
           Float64( (packet >> 16 & 0xf)) * 1.5625E-9
    
    TOT       = (packet >> 20 & 0x3ff) * 25

    dcol = (packet & 0x0FE0000000000000) >> 52                                                                  
    spix = (packet & 0x001F800000000000) >> 45                                                                    
    pix =  (packet & 0x0000700000000000) >> 44
    
    # Stefan used code unmodified. The code is for quad chip
    # but we have a single chip. Add 260 for comparison to previous runs  
    # x = Int32(dcol + pix >> 2) + 260;  # was /4                

    x = Int32(dcol + pix >> 2)                 
    y = Int32(spix + (pix & 0x3))
    
    return TOA, TOT, x, y

end # function process_chip_data_packet





"""
process global timestamp packet
returns "timemaster" whatever that is
NOT TESTED
"""
function process_global_timestamp_packet(packet)
    #------------------------------------------------------------------
    # global timestamp packet
    #------------------------------------------------------------------
    if (((packet >> 56) & 0xF) == 0x4) 
        global Timer_LSB32 = (packet >> 16) & 0xFFFFFFFF
    
    elseif (((packet >> 56) & 0xF) == 0x5)
        global Timer_MSB16 = (packet >> 16) & 0xFFFF
        # unsigned long long int timemaster
        timemaster = UInt128(Timer_MSB16)
        timemaster = (timemaster << 32) & 0xFFFF00000000
        timemaster = timemaster | Timer_LSB32;
        diff = Integer((spidrTime >> 14) - ((Timer_LSB32 >> 28) & 0x3));

        if ((spidrTime >> 14) != ((Timer_LSB32 >> 28) & 0x3))
            Timer_MSB16 = Timer_MSB16 - diff;
        end
    end

    return timemaster
end # function process_global_timestamp_packet        
