using RecipesBase
using ColorTypes
import PlotUtils: cgrad
using LaTeXStrings

@userplot Trajectories

@recipe function f(h::Trajectories)
  if length(h.args) != 2
      error("`trajectories` should be given two arguments.  Got: $(typeof(h.args))")
  end
  traj_array, sys = h.args

  @series begin
    sys.bodies
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
    size --> (800,400)
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

  body = sys.bodies[1]

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

  if typeof(field) <: VectorGridData
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
  elseif typeof(field) <: ScalarGridData
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

@userplot NozzlePlot

@recipe function f(h::NozzlePlot;fields=(),gas=DefaultPerfectGas)
  nfields = length(fields) + 1

  if nfields == 1
    length(h.args) == 1 || error("`nozzleplot` should be given one argument.  Got: $(typeof(h.args))")
    noz, = h.args
  else
    length(h.args) == 4 || error("`nozzleplot` should be given four arguments.  Got: $(typeof(h.args))")
    noz, pb, p0, T0 = h.args
    nozproc = NozzleProcess(noz,pb,p0,T0,gas=gas)
  end

  layout := (nfields,1)
  size := (500,200*nfields)
  xticks := (0:0,["Throat"])

  xline = [0,0]

  fields_lower = lowercase.(fields)

  subnum = 0
  if in("pressure",fields_lower)
    subnum += 1

    p = value.(pressure(nozproc),KPa)
    pmax = maximum(p)+50
    @series begin
      subplot := subnum
      linestyle --> :dash
      linecolor --> :black
      xline, [0,pmax]
    end
    @series begin
      subplot := subnum
      ylims --> (0,pmax)
      xlims := (-Inf,Inf)
      yguide := "Pressure (KPa)"
      #annotations := (positions(noz)[end],p[end],nozzle_quality(noz,pb,p0))
      positions(noz), p
    end
  end

  if in("temperature",fields_lower)
    subnum += 1

    T = value.(temperature(nozproc),K)
    Tmax = maximum(T)+100.0
    @series begin
      subplot := subnum
      linestyle --> :dash
      linecolor --> :black
      xline, [0,Tmax]
    end
    @series begin
      subplot := subnum

      ylims --> (0,Tmax)
      xlims := (-Inf,Inf)
      yguide := "Temperature (K)"
      positions(noz), T
    end
  end

  if in("density",fields_lower)
    subnum += 1

    ρ = value.(density(nozproc))
    ρmax = maximum(ρ)+1.0
    @series begin
      subplot := subnum
      linestyle --> :dash
      linecolor --> :black
      xline, [0,ρmax]
    end
    @series begin
      subplot := subnum

      ylims --> (0,ρmax)
      xlims := (-Inf,Inf)
      yguide := "Density (kg/m3)"
      positions(noz), ρ
    end
  end

  if in("machnumber",fields_lower) || in("mach",fields_lower)
    subnum += 1

    M = value.(machnumber(nozproc))
    Mmax = maximum(M)+1.0
    @series begin
      subplot := subnum
      linestyle --> :dash
      linecolor --> :black
      xline, [0,Mmax]
    end
    @series begin
      subplot := subnum

      ylims --> (0,Mmax)
      xlims := (-Inf,Inf)
      yguide := "Mach Number"
      positions(noz), M
    end
  end

  subnum += 1
  A = value.(areas(noz),SqCM)
  Amax = maximum(A)+10
  @series begin
    subplot := subnum
    linestyle --> :dash
    linecolor --> :black
    xline, [0,Amax]
  end
  @series begin
     subplot := subnum

     ylims --> (0,Amax)
     xlims := (-Inf,Inf)
     yguide := "Area (sq cm)"
     positions(noz),A
  end

end
