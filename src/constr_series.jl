#Se definen las funciones
function gaps_notas!(Voz::Array{Float64,2}, q::Float64) #La funcion toma de entrada el arreglo de notas y corrige los pequenios gaps que hay entre notas que terminan y que empiezan
  for i = 1:(length(Voz[:,1])-1)
      m = abs(Voz[i,2] - Voz[i+1,1])
      if m < q
          Voz[i,2] = Voz[i+1,1]
      end
  end
end

function gaps_silencios!(Voz::Array{Float64,2}) # La funcion toma de entrada el arreglo de notas y corrige los pequenios gaps que puede haber entre los silencios y donde empiezan las notas
    sm = minimum(Voz[:,2]-Voz[:,1])
    for i = 1:(length(Voz[:,1])-1) #corrige los gaps que hay entre silencios y notas
        m = mod(Voz[i,2] - Voz[i+1,1],sm)
        if m > 0
            Voz[i,2] = Voz[i,2] + (sm - m)
        end
    end
    x = length(Voz[:,1])
    n = mod(Voz[x,2] - Voz[x,1], sm) #aqui se corrige la ultima nota del arrelgo
    if n > 0
        Voz[x,2] = Voz[x,2] + (sm - n)
    end
end

function min_voces(Voces::Array{Array{Float64,2},1}, ns::Int64) #la funcion regresa el minimo valor de duracion
    mins = Array(Float64,ns)
    for i = 1:ns
        mins[i] = minimum(Voces[i][:,2] - Voces[i][:,1])
    end
    return minimum(mins)
end

function max_tempo(Voces::Array{Array{Float64,2},1}, ns::Int64) #la funcion encuentra el tiempo maximo de termino de notas, es decir donde termina la pieza
    Tfinal = Array(Float64,ns)
    for i = 1:ns
        Tfinal[i] = maximum(Voces[i][:,2])
    end
    return maximum(Tfinal)
end

function serie_notas(Voz::Array{Float64,2}, tmax::Int64) #construye la serie de tiempo teniendo como entrada el arreglo de inicio - final - voz
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
#La siguiente funcion es para filtrar csv,
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
    ind = Array(Float64,tmax)
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
function gaps_ms(s::Array{Float64,2})
    

end
