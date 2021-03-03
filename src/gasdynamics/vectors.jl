for f in (:SupersonicMachNumber,:SubsonicMachNumber,
          :SupersonicPOverP0,:SubsonicPOverP0)
  @eval $f(array::AbstractArray{T},args...;kwargs...) where T = map(x -> $f(x,args...;kwargs...),array)
end
