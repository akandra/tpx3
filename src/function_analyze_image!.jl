function analyze_image!(images, p, out1, out2, out3, out4, out5, out6, n_out)

    if p.create_image
        for i in 1 : n_out
            if out1[i] > p.mode3_t0 && out1[i] < p.mode3_tmax
                which_im = find_TOF_interval!(out3[i], p)

                if which_im > -1
                    if out6[i] > 1 && out5[i] > 1
                        images[which_im][out6[i], out5[i]] += 1
                    end
                end
            end
        end
    end
    
end

# With file 004_RealTimeTitrations_Ptlower_160C_50HzNH3_10minDosingAt100Hz_000000 from 09.02.2023, I get the folowing error:
# BoundsError: attempt to access 256×256 Matrix{Int32} at index [0, 89]
# I fix it for now by not considering these events