function get_chunk(file, data::get_word_state)
        
    if data.chunk_number == 0
        println("getting initial chunk")
        data.buffer  = reinterpret(UInt64, read(file, data.chunk*8))
    else
        println("getting chunk number ", data.chunk_number + 1)
        data.buffer  = reinterpret(UInt64, read(file, data.chunk*8))
        data.next = 1
    end

    data.chunk_number += 1
end


function write_chunk(n_out, outfile_txt, 
                    ofb1, ofb2, ofb3, ofb4, ofb5, ofb6,
                    out1, out2, out3, out4, out5, out6; 
                    write_txt=write_txt, write_bin = write_bin)
    
    if write_txt
        println("writing text file...")
        writedlm(outfile_txt, [ out1[1:n_out] out2[1:n_out] out3[1:n_out] out4[1:n_out] out5[1:n_out] out6[1:n_out] ])
    end

    if write_bin
        println("writing binary file...")
        write(ofb1, out1[1:n_out])
        write(ofb2, out2[1:n_out])
        write(ofb3, out3[1:n_out])
        write(ofb4, out4[1:n_out])
        write(ofb5, out5[1:n_out])
        write(ofb6, out6[1:n_out])
    end

end
