#Se definen las funciones
function gaps_notas!(s, q) #La funcion toma de entrada el arreglo de notas y corrige los pequenios gaps que hay entre notas que terminan y que empiezan
    for v in s
        for i = 1:(size(v)[1]-1)
            m = abs(v[i,2] - v[i+1,1])
            #printlrintln(Voz[i,2],'\t',Voz[i+1,1],'\t',m)
            if m < q
                v[i,2] = v[i+1,1]
                #println(Voz[i,2],'\t', Voz[i+1,1])
            end
        end
      end
end
#########################################################################################################################################################################
function gaps_silencios!(s,q) # La funcion toma de entrada el arreglo de notas y corrige los pequenios gaps que puede haber entre los silencios y donde empiezan las notas
    sm = q
    for v in s
        for i = 1:(size(v)[1]-1) #corrige los gaps que hay entre silencios y notas
            m = mod(v[i,2] - v[i+1,1],sm)
            if m > 0
                v[i,2] = v[i,2] + (sm - m)
            end
        end
        x = length(v[:,1])
        n = mod(v[x,2] - v[x,1], sm) #aqui se corrige la ultima nota del arrelgo
        if n > 0
            v[x,2] = v[x,2] + (sm - n)
        end
    end
end
#################################################################################################################################################
function min_voces(Voces, div) #la funcion regresa el minimo valor de duracion
    ns = size(Voces)[1]
    mins = Array{Float64}(undef,ns)
    for i = 1:ns
        dif = filter(x->x >= div,Voces[i][:,2] - Voces[i][:,1])
        mins[i] = minimum(dif[find(dif)])
    end
    return minimum(mins)
end
################################################################################################################################################
function max_tempo(Voces, ns) #la funcion encuentra el tiempo maximo de termino de notas, es decir donde termina la pieza
    Tfinal = Array{Float64}(undef,ns)
    for i = 1:ns
        Tfinal[i] = maximum(Voces[i][:,2])
    end
    return maximum(Tfinal)
end
################################################################################################################################
#La siguiente funcion corrige los numeros fraccionarios que existen en las Voces
function rounding!(voces, q)
    nv = size(voces)[1]
    for i =1:nv
        voces[i][:,1] = map(x-> ceil(x/q), voces[i][:,1])
        voces[i][:,2] = map(x-> ceil(x/q), voces[i][:,2])
    end
end
#######################################################################################################################################
function serie_notas(Voz, tmax) #construye la serie de tiempo teniendo como entrada el arreglo de inicio - final - voz
    serie = zeros(tmax)
    ini = Voz[:,1]
    dura = Voz[:,2] - Voz[:,1] #este es el arreglo de duraciones
    for i = 1:(length(dura)) #aqui se asignan las notas
        j = 0
        x = convert(Int64, ini[i]) + 1
        while dura[i] > 0
            serie[x+j] = Voz[i,3]
            j += 1
            dura[i] -= 1
        end
    end
    return(serie)
end
############################################################################################
#La siguiente funcion construye un arreglo de dos dimensiones que contiene la serie de "ritmo" y serie de notas
function series_desdoble(s::Array{Float64,2}, tmax)
    sr = Float64[]
    sp = Float64[]
    si = Float64[]
    if s[1,1] != 0 #se asegura de llenar los silencios del principio
        push!(sr,s[1,1]); push!(sp, 0); push!(si, 0)
    else
        push!(sr, s[1,2] - s[1,1]); push!(sp, s[1,3]); push!(si, s[1,4])
    end
    for i=2:(length(s[:,1])-1) #se hacen todas las notas intermedias
        if s[i,2] != s[i+1]
            push!(sr, s[i+1,1] - s[i,2]); push!(sp, 0); push!(si, 0)
        else
            push!(sr, s[i,2] - s[i,1]); push!(sp, s[i,3]); push!(si, s[i,4])
        end
    end
    push!(sr, s[end,2] - s[end,1]); push!(sp, s[end,3]); push!(si, s[end,4]) #esta es la ultima nota
    if s[end,2] < tmax
        push!(sr, tmax - s[end,2]); push!(sp, 0); push!(si, 0)  #si es que hay silencios, esto lo completa
    end
    return [sr sp si] #regresa el arreglo por ritmo, nota e intensidad.
end
####################################################################################################
#esta funcion es solo para dividir las intensidades en 0,p,mp,mf,f p.ej.
function nota_intens!(s::Array{Float64,1})
    for i=1:length(s)
        if s[i] > 0 && s[i] <= 32; s[i] = 1;
        elseif s[i] > 32 && s[i] <= 64; s[i] = 2;
        elseif s[i] > 64 && s[i] <= 96; s[i] = 3;
        elseif s[i] > 96; s[i] = 4; end
    end
end
#########################################################################################################
#La siguiente funcion es para filtrar csv, regresa una serie
#donde da el tiempo en milisegundos donde inicia la nota, termina, el pitch y la intensidad en una escala arbitraria
function filt_vozcsv(v::Array{Any,2})
    sti = Float64[]
    stf = Float64[]
    sp = Float64[]
    si = Float64[]
    for i=1:size(v)[1] #
        if v[i,6] == 0 || v[i,3] == " Note_off_c"; continue; end #aqui se salta todas las notas que "terminan"
        c = 1
        din = 0
        if v[i,5] != v[i+c,5] #&& v[i+c,6] != 0 #aqui se encuentra en donde termina la nota que inicio
            c += 1
        else
            push!(sti, v[i,2]); push!(stf, v[i+c,2]); push!(sp, v[i,5])
            if v[i,6] > 0 && v[i,6] <= 32; din = 1;
            elseif v[i,6] > 32 && v[i,6] <= 64; din = 2;
            elseif v[i,6] > 64 && v[i,6] <= 96; din = 3;
            elseif v[i,6] > 96; din = 4; end
            push!(si, din)
        end
    end
    return [sti stf sp si]
end
####################################################################################################################
function indice(tmax::Int64) #esta funcion solo regresa el arreglo que lleva el indice(numero) de las notas (o el eje x)
    ind = Array{Float64}(undef,tmax)
    for i = 1:tmax
        ind[i] = convert(Float64, i)
    end
    return(ind)
end
##############################################################################################################################
function notas_hertz!(notas::Array{Float64,1}) #esta funcion convierte el arreglo de entrada de numero de nota a su valor en hertz con A4 = 440Hz
    a4 = 440.0
    x = 1.0/12.0
    r = 2 ^ x
    for i = 1:length(notas)
        if notas[i] == 0.0 #si la nota ya es silencio, la omite
            continue
        end
        d = notas[i] - 69 # 69 es el A4
        notas[i] = a4 * r ^ d
    end
end
#################################################################################
function note_vec(v::Array{Float64,2})
    n = size(v)[1]
    vec = AbstractString[]
    for i=1:n
        nota = string(v[i,1], ',', v[i,2], ',', v[i,3])
        push!(vec, nota)
    end
    return vec
end
#######################################################################################
function csvtoserie(s, sd)
    #sd is the subdivision, is the smallest time duration of a note, usualy is 8 (1/8) or 16 (1/16), you can put 8 to try.
    mq = s[1,6] #quarter of a note
    nv = s[findlast(s[:,3], " Note_on_c"),1] - 1 #estimates how many voices are in the midi
    voces = Array{Matrix}(undef,nv) #initialze an array
    if findfirst(s[:,3], " Note_off_c") == 0 #checks if the midi has events of note_off
        for i = 2:(nv+1) #if does not, it construct the series in this way
            b = s[s[:,1].==i,:] #takes the events of the voice i
            if size(b[b[:,3].==" Note_on_c", :])[1] == 0; continue; end #checks if there are notes in the channel
            voces[i-1] = float(filt_vozcsv(b[b[:,3].==" Note_on_c", :])) #construct an array of information of initial time, finish time, pitch and intensity.
        end
    else
        for i = 2:(nv+1) #if has note_off events , it does this way
            ini = findfirst(s[s[:,1].==i,3], " Note_on_c") #initial time
            fin = findlast(s[s[:,1].==i,3], " Note_off_c") #finish time
            if ini == 0 || fin == 0; continue; end
            voces[i-1] = float(filt_vozcsv(s[s[:,1].==i,:][ini:fin,:])) #construct an array of information of initial time, finish time, pitch and intensity
        end
    end
    #voces = voces[map(x -> isdefined(voces, x ), 1:length(voces))] #these two steps are just to trow away the voices that are empty.
    filter!(x -> size(x)[1] != 0,voces)
    nv = size(voces)[1] #the real number of voices
    q = mq/ round(Int, mq / min_voces(voces, mq / sd)) #defines the unit of time  q
    gaps_notas!(voces,q) #fixes some gaps of miliseconds between the notes
    gaps_silencios!(voces,q) #fixes some gaps of milisecons of rests
    rounding!(voces, q) #rounds up some possible decimals in the time
    tmax = round(Int,max_tempo(voces,nv)) #gets the total lenght of the piece
    series = zeros(tmax,nv+1)  #this is the initialization of the output array
    series[:,1] = indice(tmax) #this is just the index
    for i=2:(nv+1) #constructs the output time series
        series[:,i] = serie_notas(voces[i-1],tmax)
    end
    return series
end
##############################################################################
function csvtointervs(s::Array{Any,2}, sd::Int64)
    #sd is the subdivision, had to correct it
    mq = s[1,6] #quarter of a note
    nv = s[findlast(s[:,3], " Note_on_c"),1] - 1
    voces = Array(Array{Float64,2},nv)
    if findfirst(s[:,3], " Note_off_c") == 0
        for i = 2:(nv+1)
            b = s[s[:,1].==i,:]
            #ini = findfirst(s[s[:,1].==i,3], " Note_on_c")
            #fin = findlast(s[s[:,1].==i,3], " Note_on_c")
            #println(ini,'\t', fin)
            if size(b[b[:,3].==" Note_on_c", :])[1] == 0; continue; end
            voces[i-1] = float(filt_vozcsv(b[b[:,3].==" Note_on_c", :]))
        end
    else
        for i = 2:(nv+1)
            ini = findfirst(s[s[:,1].==i,3], " Note_on_c")
            fin = findlast(s[s[:,1].==i,3], " Note_off_c")
            #println(ini,'\t', fin)
            if ini == 0 || fin == 0; continue; end
            voces[i-1] = float(filt_vozcsv(s[s[:,1].==i,:][ini:fin,:]))
        end
    end
    voces = voces[map(x -> isdefined(voces, x ), 1:length(voces))]
    filter!(x -> size(x)[1] != 0,voces)
    nv = size(voces)[1]
    q = mq/ round(Int, mq / min_voces(voces, mq / sd))
    gaps_notas!(voces,q)
    gaps_silencios!(voces,q)
    rounding!(voces, q)
    inters = Array(Array{Any,1},nv)
    for i=1:(nv)
        inters[i] = inter_symbol(voces[i])
    end
    return inters
end
#############################################################################
function inter_symbol(s::Array{Float64,2})
    nuevas = Array{Float64,1}()
    if s[1,1] != 0.0; push!(nuevas, 1000.0); end
    push!(nuevas,s[1,3])
    for i = 1:(size(s)[1]-1)
        if s[i,2] != s[i+1,1]
            push!(nuevas,1000.0)
        else
            push!(nuevas,s[i+1,3])
        end
    end
    inters = Array{Float64,1}()
    for i = 1:(length(nuevas)-1)
        if nuevas[i] > 500; push!(inters, nuevas[i])
        elseif nuevas[i+1] > 500; continue;
        else
            push!(inters, nuevas[i+1]-nuevas[i])
        end
    end
    return convert(Array{Any,1},inters)
end
###############################################################################
function rank_series(rf, s)
    inte=convert(Array{Float64,1},s)
    rs = zeros(length(s))
    for i = 1:size(rf)[1]
        indx = findin(inte,rf[i,2])
        for j = 1:length(indx)
            rs[indx[j]] = i
        end
    end
    return rs
end

function rank_series_i(rf, s)
    inte=convert(Array{Float64,1},s)
    rs = zeros(length(s))
    N = size(rf)[1]
    for i = 1:size(rf)[1]
        #if rf[i,2] == 0; continue; end
        indx = findin(inte,rf[i,2])
        for j = 1:length(indx)
            rs[indx[j]] = -i + N + 1
        end
    end
    return rs
end

function freq_series(rf, s)
    inte=convert(Array{Float64,1},s)
    rs = zeros(length(s))
    for i = 1:size(rf)[1]
        indx = findin(inte,rf[i,2])
        for j = 1:length(indx)
            rs[indx[j]] = rf[i,3]
        end
    end
    return rs
end
################################################################################
function rank_seriestxt(rf::Array{Any,2}, s::Array{AbstractString,1}, lt::Array{AbstractString,1})
    rs = zeros(length(s))
    for i = 1:length(s)
        for k =1:length(lt)
            if contains(s[i],lt[k])
                rs[i] = k; break
            end
        end
    end
    return rs
end
################################################################################
function red_hv(rf::Array{Any,2},rs::Array{Float64,1},am::Array{Float64,2})
    n = size(rf)[1]
    rsn = convert(Array{Int64,1}, rs)
    red_m = zeros(n,n)
    for i = 1:size(am)[1]
        for j = 1:size(am)[2]
            if am[i,j] == 1; red_m[rsn[i],rsn[j]] = 1; end
        end
    end
    return red_m
end
###############################################################################
function consonance_series(s::Array{Float64,2})
    nv = size(s)[2]
    n = size(s)[1]
    CR = [1,11,8,6,4,3,12,2,7,5,9,10] #este es el rango de consonancia
    cs = zeros(n)
    for k = 1:n
        c = 0
        for i = 2:(nv-1)
            for j = (i+1):nv
                if s[k,i] == 0 || s[k,j] == 0; continue; end #aqui asegura de dejar en cero la consonancia entre una nota y nada mas
                cn = CR[Int(mod(abs(s[k,i] - s[k,j]),12) + 1)]
                if cn > c; c = cn; end
            end
        end
        cs[k] = c
    end
    return [s[:,1] cs]
end
################################################################################
function note_serie(s::Array{Any,2})
    nv = s[findlast(s[:,3], " Note_on_c"),1] - 1
    voces = Array(Array{Float64,2},nv)
    if findfirst(s[:,3], " Note_off_c") == 0
        for i = 2:(nv+1)
            b = s[s[:,1].==i,:]
            #ini = findfirst(s[s[:,1].==i,3], " Note_on_c")
            #fin = findlast(s[s[:,1].==i,3], " Note_on_c")
            #println(ini,'\t', fin)
            if size(b[b[:,3].==" Note_on_c", :])[1] == 0; continue; end
            voces[i-1] = float(filt_vozcsv(b[b[:,3].==" Note_on_c", :]))
        end
    else
        for i = 2:(nv+1)
            ini = findfirst(s[s[:,1].==i,3], " Note_on_c")
            fin = findlast(s[s[:,1].==i,3], " Note_off_c")
            #println(ini,'\t', fin)
            if ini == 0 || fin == 0; continue; end
            voces[i-1] = float(filt_vozcsv(s[s[:,1].==i,:][ini:fin,:]))
        end
    end
    return voces
end
################################################################################
function interval_serie(s::Array{Array{Float64,2},1})
    nv = size(s)[1]
    interval = Array(Array{Any,1},nv)
    for i = 1:nv
        interval[i] = diff(s[i][:,3])
    end
    return interval
end
function interval_serie(s::Array{Array{Float64,2},1},n::Int64)
    nv = size(s)[1]
    interval = Array(Array{Any,1},nv)
    for i = 1:nv
        if n == 1
            interval[i] = diff(shuffle(s[i][:,3]))
        end
    end
    return interval
end
################################################################################
#Next function will return the matrix of horizontal visibility graph
function h_visibility(s)
    adj_mat = zeros(length(s),length(s))
    for i=1:(length(s)-2)
        adj_mat[i,i+1] = 1; adj_mat[i+1,i] = 1
        for j=(i+2):length(s)
            kt = maximum(s[(i+1):(j-1)])
            if s[i] > kt && s[j] > kt
                adj_mat[i,j] = 1; adj_mat[j,i] = 1
                if s[j] >= s[i]; break; end
            end
        end
    end
    adj_mat[end-1,end] = 1; adj_mat[end,end-1] = 1
    return adj_mat
end


function visibility(s)
    adj_mat = zeros(length(s),length(s))
    for i=1:(length(s)-2)
        adj_mat[i,i+1] = 1; adj_mat[i+1,i] = 1
        for j=(i+2):length(s)
            k = findlast(adj_mat[i,:]) #locates the index with the last link to i
            m1 = (s[k]-s[i])/(k-i); m2 = (s[j]-s[k])/(j-k) #estimate the slopes
            if s[k] <= s[j] &&  m1 < m2  #using the criteria to link the nodes
                adj_mat[i,j] = 1; adj_mat[j,i] = 1
            end
        end
    end
    adj_mat[end-1,end] = 1; adj_mat[end,end-1] = 1
    return adj_mat
end

function hd_visibility(s::Array{Int64,1})
    adj_mat = zeros(length(s),length(s))
    for i=1:(length(s)-2)
        adj_mat[i,i+1] = 1; adj_mat[i+1,i] = 1
        for j=(i+2):length(s)
            if s[i] > maximum(s[(i+1):(j-1)]) && s[j] > maximum(s[(i+1):(j-1)])
                adj_mat[i,j] = 1;
            end
        end
    end
    adj_mat[end-1,end] = 1; adj_mat[end,end-1] = 1
    return adj_mat
end

function hd_visibility(s::Array{Float64,1})
    adj_mat = zeros(length(s),length(s))
    for i=1:(length(s)-2)
        adj_mat[i,i+1] = 1; adj_mat[i+1,i] = 1
        for j=(i+2):length(s)
            if s[i] > maximum(s[(i+1):(j-1)]) && s[j] > maximum(s[(i+1):(j-1)])
                adj_mat[i,j] = 1;
            end
        end
    end
    adj_mat[end-1,end] = 1; adj_mat[end,end-1] = 1
    return adj_mat
end
################################################################################
#next function constructs a series of blocks from a time series, by the criteria
#of the HVG
function hv_blocks(s)
    out = Array{Array{Float64,1},1}()
    for i=1:(length(s)-2)
        b = 0
        for j=(i+2):length(s)
            kt = maximum(s[(i+1):(j-1)])
            if s[i] > kt && s[j] > kt
                b = j
            end
        end
        if b != 0
            push!(out,s[i:b])
        else
            b == 0; push!(out,s[i:i+1])
        end
    end
    return out
end
################################################################################
function v_blocks(s)
    adj_mat = zeros(length(s),length(s))
    for i=1:(length(s)-2)
        adj_mat[i,i+1] = 1; adj_mat[i+1,i] = 1
        for j=(i+2):length(s)
            k = findlast(adj_mat[i,:]) #locates the index with the last link to i
            m1 = (s[k]-s[i])/(k-i); m2 = (s[j]-s[k])/(j-k) #estimate the slopes
            if s[k] <= s[j] &&  m1 < m2  #using the criteria to link the nodes
                adj_mat[i,j] = 1; adj_mat[j,i] = 1
            end
        end
    end
    adj_mat[end-1,end] = 1; adj_mat[end,end-1] = 1
    out = Array{Array{Float64,1},1}()
    for i = 1:(length(s)-1)
        k = findlast(adj_mat[i,:])
        push!(out, s[i:k])
    end
    return out
end
function vr_blocks(s) #recursive blocks
    adj_mat = zeros(length(s),length(s))
    for i=1:(length(s)-2)
        adj_mat[i,i+1] = 1; adj_mat[i+1,i] = 1
        for j=(i+2):length(s)
            k = findlast(adj_mat[i,:]) #locates the index with the last link to i
            m1 = (s[k]-s[i])/(k-i); m2 = (s[j]-s[k])/(j-k) #estimate the slopes
            if s[k] <= s[j] &&  m1 < m2  #using the criteria to link the nodes
                adj_mat[i,j] = 1; adj_mat[j,i] = 1
            end
        end
    end
    adj_mat[end-1,end] = 1; adj_mat[end,end-1] = 1
    out = Array{Array{Float64,1},1}()
    for i = 1:(length(s)-1)
        k = findfirst(adj_mat[i,:])
        push!(out,s[i:k])
        while k < findlast(adj_mat[i,:])
            k = findnext(adj_mat[i,:],k+1)
            if isempty(s[i:k]); continue; end
            push!(out, s[i:k])
        end
    end
    return out
end
################################################################################
#returns a list of edges and a list of nodes.
function admat_out(am)
    mat_out = []; n_out = []; ng = size(am)[1]
    for i=1:ng #this step is for printing the output transition matrix.
        for j=1:ng
            if am[i,j]==0; continue; end
            push!(mat_out, [i j])
            push!(n_out, i)
        end
    end
    #writecsv("smw-a1.csv", vcat(mat_out...))
    #writecsv("smw-a1n.csv", vcat(n_out...))
    return n_out, mat_out
end
################################################################################
function vmin_blocks(s) #minimum size blocks. with VG algorithm
    out = Array{Array{Float64,1},1}()
    n = 1
    while n < length(s)
        br = false
        if n==(length(s)-1)
            b = Array{Float64,1}()
            push!(b,s[n]); push!(b,s[end])
            push!(out,b)
            n = length(s)
        end
        for i=n:(length(s)-2)
            #println("i=",i)
            b = Array{Float64,1}()
            push!(b, s[i]);push!(b,s[i+1])
            k=i+1
            for j=(i+2):length(s)
                val = s[j] + (s[i]-s[j])*(j-k)/(j-i)
                #println("val =",val,'\t',"s[k]=",s[k])
                if s[k] < val
                    push!(b, s[j])
                    k+=1
                else
                    n = j-1
                    #println("n=",n)
                    br = true
                    break
                end
                if j==length(s); br=true;n=length(s);end
            end
            if br
                #println(b)
                push!(out, b)
                break
            end

        end
    end
    return out
end
################################################################################
function hvi_blocks(s::Array{Float64,1})
    out = Array{Array{Int64,1},1}()
    for i=1:(length(s)-2)
        b = 0
        for j=(i+2):length(s)
            if s[i] > maximum(s[(i+1):(j-1)]) && s[j] > maximum(s[(i+1):(j-1)])
                b = j
            end
        end
        if b != 0
            push!(out,collect(i:b))
        else
            b == 0; push!(out,[i,i+1])
        end
    end
    return out
end
function hv_links(s)
    adj_mat = zeros(length(s),length(s))
    for i=1:(length(s)-2)
        adj_mat[i,i+1] = 1; adj_mat[i+1,i] = 1
        for j=(i+2):length(s)
            kt = maximum(s[(i+1):(j-1)])
            if s[i] > kt && s[j] > kt
                for k = (i+1):j
                    adj_mat[i,k] = 1; adj_mat[k,i] = 1
                end
                break
            end
                #break
#            elseif j == length(s)
#                for k = (i+1):j
#                    adj_mat[i,k] = 1; adj_mat[k,i] = 1
#                end
            #end
        end

    end
    adj_mat[end-1,end] = 1; adj_mat[end,end-1] = 1
    return adj_mat
end

#################################################################################
#Next function classifies blocks with the hamming distances and k-medoids method
function group_blocks(b::Array{Array{Float64,1},1}, pen::Float64)
    mx = maximum(map(length, b))
    groups = Array{Array{Array{Float64,1},1},1}()
    for l = 2:mx
        bl = filter(x -> length(x) == l, b)
        if isempty(bl); continue; end #checking for an empty array.
        n = length(bl)
        mat_d = zeros(n,n)
        for i= 1:(n-1)
            for j = (i+1):(n)
                mat_d[i,j] = mat_d[j,i] = euclidean(bl[i],bl[j]) #estimating euclidean distances
            end
        end
        err = 10
        k = 1
        out = kmedoids(mat_d,k)
        while err > pen
            out = kmedoids(mat_d, k)
            err = maximum(out.acosts)
            k += 1
        end
        nm = length(out.medoids)
        meds = out.medoids
        ass = out.assignments
        for i=1:nm
            group = Array{Array{Float64},1}()
            ind_ass = find(x->x==i, ass)
            for j=1:length(ind_ass)
                push!(group, bl[ind_ass[j]])
            end
            push!(groups,group)
        end
    end
    return groups

end
################################################################################
function hgroup_blocks(b::Array{Array{Float64,1},1}, pen::Float64)
    n = size(b)[1]
    mat_d = zeros(n,n)
    for i= 1:(n-1)
        for j = (i+1):(n)
            mat_d[i,j] = mat_d[j,i] = ham_distexp(b[i],b[j])
        end
    end
    err = 4
    k = 1
    out = kmedoids(mat_d,k)
    while err > pen
        out = kmedoids(mat_d, k)
        err = maximum(out.acosts)
        k += 1
    end
    nm = length(out.medoids)
    groups = Array{Array{Array{Float64,1},1},1}()
    meds = out.medoids
    ass = out.assignments
    for i=1:nm
        group = Array{Array{Float64},1}()
        ind_ass = find(x->x==i, ass)
        for j=1:length(ind_ass)
            push!(group, b[ind_ass[j]])
        end
        push!(groups,group)
    end
    return groups
end
################################################################################
#next function returns the original values of the timeseries, given blocks of rank frequency.
function rf_blocks(b::Array{Array{Float64,1},1}, rf::Array{Any,2})
    bout = Array{Array{Float64,1},1}()
    for i = 1:length(b)
        b1 = Array{Float64,1}()
        for j = 1:length(b[i])
            push!(b1, rf[Int64(b[i][j]),2])
        end
        push!(bout, b1)
    end
    return bout
end
##################################################################################
#next function inverts the time series
function inv_serie(s)
    s_in = map(x-> -x + maximum(s) + minimum(s), s)
    return s_in
end
################################################################################
#this function returns a symbolic series from text file.
function txt_series(f)
    #read the file
    s = read(f, String)
    #st2 = replace(s2,r"\n"i,"")
    #st2 = replace(st2,r"\r"i,"")
    #next lines are for removing every character that is not a letter.
    st = replace(s, r"[^a-z Ñ ñ Á á Ĉ ĉ Ĝ ĝ Ĥ ĥ É é Í í Ĵ ĵ Ŝ ŝ Ŭ ŭ Ó ó Ú ú ü Ä ä Ü ü Ö ö ẞ ß]"i => " ")
    st = replace(st, r"[A-Z Ñ Á É Í Ó Ú Ĉ Ĝ Ĥ Ĵ Ŝ Ŭ Ä Ü Ö ẞ]", lowercase) #replace(s, r"[^A-Z]" => lowercase)
    #st = replace(s, r"[^a-z]"i => " ")
    #st = replace(st, r"[A-Z]", lowercase) #replace(s, r"[^A-Z]" => lowercase)

    l = convert(Array{AbstractString,1},split(st," "))
    filter!(x -> x!="",l)
    s = split(join(l," "),"")
    filter!(isascii, s)
    return s
end
