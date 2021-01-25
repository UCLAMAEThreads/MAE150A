### Plotting arrows on streamlines

"""
    get_segment_coords(p::Plot,level::Integer)

Get the x,y coordinates of contour level `level` in plot `p`.
"""

function get_segment_coords(ps::Plots.Plot{Plots.PyPlotBackend},level::Integer)
    contourset = get_series_handle(ps)

    if level > length(contourset."collections")
      error("This contour level is not on the list.")
    end
    lc = get(contourset."collections",level-1)
    seg_array_list = PyCall.pycall(lc."get_segments",Array)

    return seg_array_list
end

function get_series_handle(ps::Plots.Plot{Plots.PyPlotBackend})
    sl = ps.series_list[1]
    contourset = sl.plotattributes[:serieshandle][1]
    return contourset
end

function arrowhead_coords(x,y,u,v,scale;open_angle::Float64=π/3)
    R = [cos(0.5*open_angle) sin(0.5*open_angle);
          -sin(0.5*open_angle) cos(0.5*open_angle)]

    uvec = [u,v]/sqrt(u^2+v^2)

    xhead = [x,y]
    x1 = xhead .- scale*(R*uvec)
    x2 = xhead .- scale*(R'*uvec)

    xdata = [x1[1],x,x2[1]]
    ydata = [x1[2],y,x2[2]]
    return xdata, ydata
end

function arrowhead!(ps::Plots.Plot{Plots.PyPlotBackend},x,y,u,v,scale;open_angle::Float64=π/3)

    xdata, ydata = arrowhead_coords(x,y,u,v,scale,open_angle=open_angle)

    return Plots.plot!(ps,xdata,ydata,color=:black,linewidth=1)
end

function add_arrow!(ps::Plots.Plot{Plots.PyPlotBackend},seg::Array{Float64,2},elements;num_arrows=1)

    if size(seg,1) <= 1
        return ps
    end

    lastrow = size(seg,1)
    if seg[end,:] == seg[end-1,:]
        lastrow -= 1
    end
    firstrow = 1
    if lastrow-firstrow >= 2
      firstrow += 1
    elseif lastrow-firstrow <= 0
      return ps
    end
    deg = (lastrow-firstrow) < 3 ? 1 : 3

    spl = ParametricSpline(seg[firstrow:lastrow,:]',k=deg)

    interv = 1/(num_arrows+1)
    for f in interv:interv:1-interv
        xhead, yhead = spl(f)
        w = PotentialFlow.induce_velocity(xhead+im*yhead,elements,0)
        arrowhead!(ps,xhead,yhead,real(w),imag(w),0.1)
    end
    return ps
end

"""
    add_arrow!(p::Plot,level::Integer,elements,[num_arrows=2])

Add arrow(s) on a given contour level of the plot `p`. The `elements` are the potential flow element(s)
that generate the flow. The number of arrows can be specified with the optional argument, but defaults
to 2
"""
function add_arrow!(ps::Plots.Plot{Plots.PyPlotBackend},level::Integer,el;a...)
    seg_list = get_segment_coords(ps,level)
    if length(seg_list) > 0
        for seg in seg_list
            add_arrow!(ps,seg,el;a...)
        end
    else
        return ps
    end
end

function add_arrows!(ps::Plots.Plot{Plots.PyPlotBackend},el;a...)
    contourset = get_series_handle(ps::Plots.Plot)
    nlev = length(contourset."levels")
    for lev in 1:2:nlev
        add_arrow!(ps,lev,el;a...)
    end
    return ps
end
