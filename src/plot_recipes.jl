using RecipesBase
using ColorTypes
import PlotUtils: cgrad, palette, color_list
using LaTeXStrings

@userplot Trajectories

@recipe function f(h::Trajectories)
  if length(h.args) != 2
      error("`trajectories` should be given two arguments.  Got: $(typeof(h.args))")
  end
  traj_array, sys = h.args

  @series begin
    sys.base_cache.bl
  end

  for traj in traj_array
    @series begin
      linewidth := 2
      traj[1,:], traj[2,:]
    end
  end

  @series begin
    xlims --> (-1,3)
    ylims --> (-1,1)
    xguide := L"x"
    yguide := L"y"
    title := "Particle trajectories"
    aspect_ratio := 1
    size --> (700,400)
    ()
  end


end


@userplot FieldTrajectory

@recipe function f(h::FieldTrajectory;fieldlabel="Field",deriv=0)
  if length(h.args) != 3
      error("`fieldtrajectory` should be given three arguments.  Got: $(typeof(h.args))")
  end
  traj, field, sys = h.args

  xtraj = traj[1,:]
  ytraj = traj[2,:]

  body = sys.base_cache.bl[1]

  layout := (2,1)
  size --> (600,600)
  xlims := (-1,2)

  @series begin
    subplot := 1
    ylims --> (-0.5,0.5)
    xguide := L"x"
    yguide := L"y"
    aspect_ratio := 1
    body
  end

  @series begin
    subplot := 1
    aspect_ratio := 1
    ylims --> (-0.5,0.5)
    title := "Particle trajectory"
    xtraj, ytraj
  end

  yb = -20
  yt = 20
  Xr = [minimum(body.x), maximum(body.x), maximum(body.x), minimum(body.x), minimum(body.x)]
  Yr = [yb,yb,yt,yt,yb]

  @series begin
    subplot := 2
    fillrange := 0
    fillalpha --> 0.2
    linealpha --> 0.2
    fillcolor := :lightgray
    linecolor := :lightgray
    label := "Body location"
    Xr,Yr
  end

  if typeof(field) <: ViscousFlow.VectorGridData
    utraj,vtraj = field_along_trajectory(field,sys,traj,deriv=deriv)
    minuv = min(minimum(utraj),minimum(vtraj)) - 0.5
    maxuv = max(maximum(utraj),maximum(vtraj)) + 0.5

    @series begin
      subplot := 2
      xguide := L"x"
      label := L"$x$ %$fieldlabel"
      xtraj, utraj
    end

    @series begin
      subplot := 2
      xguide := L"x"
      label := L"$y$ %$fieldlabel"
      title := "$fieldlabel components along trajectory"
      legend := true
      ylims := (minuv,maxuv)
      xtraj, vtraj
    end
  elseif typeof(field) <: ViscousFlow.ScalarGridData
    straj = field_along_trajectory(field,sys,traj,deriv=deriv)

    @series begin
      subplot := 2
      xguide := L"x"
      label := "$fieldlabel"
      title := "$fieldlabel along trajectory"
      ylims := (minimum(straj)-0.5,maximum(straj)+0.5)
      legend := true
      xtraj, straj
    end
  end

end
