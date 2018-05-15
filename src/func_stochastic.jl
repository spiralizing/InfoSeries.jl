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
################################################################################
function get_glinks(g, m)
    ng = size(g)[1]
    indx = Array{Int64,1}()
    for i = 1:ng
        if isempty(find(x->x[1]==m,g[i]))==false
            push!(indx, i)
        end
    end
    return indx
end

function get_slinks(g, m)
    ng = size(g)[1]
    indx = Array{Int64,1}()
    for i = 1:ng
        if isempty(find(x->x[end]==m,g[i]))==false
            push!(indx, i)
        end
    end
    return indx
end
################################################################################
#Next function is for constructing a small world-like network from a set of groups of blocks.
function init_swmat(i_groups, ϵ)
    ng = size(i_groups)[1]
    ind_g = collect(1:ng)
    am = zeros(ng,ng)
    lv = Array{Int64,1}(); push!(lv, 1) #to keep track of the nodes visited
    il = 1
    nl = length(lv)
    counter = 0
    counter2 = 0
    tam_res = Array{Int64,1}()
    stop=false
    res_nodes = Array{Int64,1}()
    while stop ==false
        lg = get_glinks(i_groups[1:end .!= il],i_groups[il][1][end]) #this locates the possible transitions
        if isempty(lg)
            println("There are no possible transitions")
            #break
            println("Groups ",il, " from ", ng)
            counter2 += 1
            if counter2 > 3; stop=true; end
            continue
        end
        l = rand(lg) #chooses one randomly
        c_v = findin(lv,l) #is checking if the node was already visited
        #println("Transition from ",il," to ",l)
        if isempty(c_v) #in the case that the node has not been visited.
            am[il,l] = 1 #creates a link between the nodes il -> l
            il = l #now we will move from the new node.
            push!(lv, il) #adding the node to the visited nodes.
        else #now for the case that the node has already been visited
            while length(lg) > 0
                #println("possible transitions", '\t', lg, "\t node visited before \t",l)
                filter!(x->x!=l,lg) #removes the node that has been visited
                if isempty(lg)
                    counter += 1
                    res = setdiff(ind_g,lv)
                    push!(tam_res, length(res))
                    if counter > ng
                        #println(counter)
                        if length(unique(tam_res[end-length(res):end])) == 1
                            stop=true
                            res_nodes = res
                        end
                    end
                    #println(res)
                    il = rand(res)
                    break
                end
                l = rand(lg) #a new node is chosen
                c_v = findin(l,lv)
                if isempty(c_v)
                    am[il,l] = 1
                    il = l
                    push!(lv, il)
                    break
                #else
                #    println(lv[c_v],'\t',l,'\t',lg)
                end
            end
            #println("nodes visited:",'\t',nl)
            #break
        end
        nl = length(lv)
        #println(nl,'\t',"nodes visited")
    end
    for i = 1:length(res_nodes)
        lg = get_glinks(i_groups[1:end .!= res_nodes[i]],i_groups[res_nodes[i]][1][end])
        if isempty(lg); continue;end
        l = rand(lg)
        am[res_nodes[i],l] = 1
    end
    for i = 1:ng
        if isempty(find(x-> x==1, am[i,:]))
            #println("El nodo ",i," no tiene salida")
            salidas = get_glinks(i_groups[1:end .!=i],i_groups[i][1][end])
            if isempty(salidas); continue; end
            el = rand(get_glinks(i_groups[1:end .!=i],i_groups[i][1][end]))
            am[i,el] = 1
        end
        if isempty(find(x-> x==1, am[:,i]))
            #println("El nodo ",i," no tiene entrada")
            entradas = get_slinks(i_groups[1:end .!=i],i_groups[i][1][1])
            if isempty(entradas); continue;end
            el = rand(get_slinks(i_groups[1:end .!=i],i_groups[i][1][1]))
            am[el,i] = 1
        end
    end
    for i = 1:ng #this step is just to grow the number of edges.
        if rand() < ϵ
            pl = find(x->x==0,am[i,:])
            lg = get_glinks(i_groups[pl[1:end .!= i]],i_groups[i][1][end])
            if isempty(lg); continue;end
            link = rand(lg)
            am[i,link] = 1
        end #initialize the network, links.
    end
    return am
end
