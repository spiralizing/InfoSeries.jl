#Aqui hay funciones para procesos estocasticos.
#############################################################
function mat_trans1(s::Array{Float64,1},tam::Int64)
    M = zeros(tam,tam)
    for i=1:(length(s)-1)
        M[s[i]+1,s[i+1]+1] = 1
    end
    return M
end
function mat_trans_p(s::Array{Float64,1},tam::Int64)
    M = zeros(tam,tam)
    for i=1:(length(s)-1)
        M[s[i]+1,s[i+1]+1] += 1
    end
    for i=1:128
        if sum(M[i,:]) == 0; continue; end
        M[i,:] = M[i,:]/sum(M[i,:])
    end
    return M
end
function mat_trans(s::Array{Float64,1})
    s = convert(Array{Int64,1}, s)
    tam = maximum(s)
    M = zeros(tam,tam)
    for i=1:(length(s)-1)
        M[s[i],s[i+1]] = 1
    end
    return M
end
