function update_fields!(fields, fds, clock, frame; skip_update)

    fieldnames = Symbol.(keys(fds.fields))
    for fieldname in fieldnames
        fieldname ∈ skip_update && continue
        set!(fields[fieldname], fds[fieldname][frame])
    end
    
    # Set previous state
    f2 = max(frame - 1, 1)
    
    for fieldname in fieldnames
        nextfieldname = Symbol(fieldname, :_prev)
        nextfieldname ∈ skip_update && continue
        nextfieldname ∉ keys(fields) && continue
        set!(fields[nextfieldname], fds[fieldname][f2])
    end

    return nothing
end
