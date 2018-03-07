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
################################################################################
#Next functions are different kind of random walks over a matrix
#arguments are: the adjacency matrix or normalized transition matrix,
#Starting point, if not, leave blank, and number of steps.
function r_walk(adj_mat, st, tmax)
    s = Array{Float64,1}()
    if sum(adj_mat[:,1] > 2.0) #this is to check if the matrix is normalized or not
        for i = 1:tmax
            try
                p = find(adj_mat[:,st]) #encuentra las transiciones posibles
                nv = rand(p)
                push!(s,nv)
                st = nv
            catch
                println("Error with the adjacency matrix")
            end
        end
    else
        for i = 1:tmax
            p = find(adj_mat[st,:]) #encuentra las transiciones posibles
            pc = cumsum(prob_mat[st,p]) #hace un arreglo con la distrubucion de probabilidad acumulada en las transiciones.
            α = rand() #se tira una moneda entre 0-1
            pi = findfirst(sort([α;pc][:,1]), α) #encuentra el lugar de transicion
            #v = rand(groups[p[pi]]) #se escoge uno de los bloques dentro del grupo
            push!(mnotas,p[pi]) #se agrega a la secuencia solo la primer nota del bloque elegido
        end
    end
end


##################################################################################
#this function finds the blocks that start with r in the g array of blocks
function findgroup(g, r)
    ng = size(g)[1]
    ind = []
    for i = 1:ng
        push!(ind,isempty(find(x-> x[1]==r, g[i])))
    end
    return find(x -> x==false, ind)
end
