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
function series_desdoble(s::Array{Float64,2})
    sr = Float64[]
    sp = Float64[]
    if s[1,1] != 0
        push!(sr,s[1,1])
        push!(sp, 0)
    else
        push!(sr, s[1,2] - s[1,1])
        push!(sp, s[1,3])
    end
    for i=2:length(s[:,1])
        push!(sr, s[i,2] - s[i,1])
        push!(sp, s[i,3])
    end
    return [sr sp]
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
