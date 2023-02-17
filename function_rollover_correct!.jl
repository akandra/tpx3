"""
    correct timestamps for rollover
        arguments:
            time            timestamp time      {Float}
            prior_time      previous timestamp  {Float}
            rollover        rollover time       {Float}
            time_correction correction to add   {Float}
         modifies in place:
            time, prior_time, correction
"""
function rollover_correct!(time, prior_time, rollover, correction)
    
    time += correction
    
    if time - prior_time < -10.0
        time       += rollover
        correction += rollover
    end

    if time - prior_time > 10.0
        time       -= rollover
        correction -= rollover
    end
    
    
    # println("time, prior: ", time, " ", prior_time)   
    prior_time  = time

    return time, prior_time, correction

end # function rollover_correct!