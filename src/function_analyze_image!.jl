function analyze_image!(image::Matrix{Int32}, 
                        p::pars, 
                        out1::Vector{Float32}, 
                        out2::Vector{Float32}, 
                        out3::Vector{Float32}, 
                        out4::Vector{UInt16}, 
                        out5::Vector{UInt8}, 
                        out6::Vector{UInt8}, 
                        n_out::Int)

    for i in 1:n_out
        if out3[i] > p.tof_gate_min && out3[i] < p.tof_gate_max
            image[out6[i], out5[i]] += 1
        end
    end
    
end

function analyze_image_kt!(image, out2, out3, TOFgate, kt_gate, out5, out6, n_out)

    ta = now()

    for i in 1:n_out
        if out3[i] > TOFgate[1] && out3[i] < TOFgate[2]
        if out2[i] > kt_gate[1] && out2[i] < kt_gate[2]
            image[out6[i],out5[i]] += 1
        end
        end
    end
    
    return now() - ta
end