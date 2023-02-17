function analyze_image!(images, p, out1, out2, out3, out4, out5, out6, n_out)

    if p.create_image
        for i in 1 : n_out
            which_im = find_TOF_interval!(out3[i], p)

            if typeof(which_im) != Nothing
                images[which_im][out6[i], out5[i]] += 1
            end
            
        end
    end
    
end

