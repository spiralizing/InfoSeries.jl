#estas son las funciones que utiliza mi programa que calcula el nDFA
#####################################################################################
#Esta funcion integra la serie x en yᵢ = sum [xᵢ - <x>]
function integrate_serie!(serie::Array{Float64,1})
    n = length(serie)
    out = zeros(n)
    for i = 1:n
        s = 0
        for j = 1:i
            s += serie[j] - mean(serie)
        end
        out[i] = s
    end
    serie = out
end

function integrate(serie::Array{Float64,1})
    n = length(serie)
    out = zeros(n)
    for i = 1:n
        s = 0
        for j = 1:i
            s += serie[j] - mean(serie)
        end
        out[i] = s
    end
    return out
end
################################################################################
#esta funcion regresa la serie de magnitudes a partir de una serie original
function magnitude_serie(serie::Array{Float64,1})
    mag = zeros(Float64,length(serie))
    for i = 1:(length(serie))
        mag[i] = abs(serie[i])
    end
    return mag - mean(mag)
end
################################################################################
#esta funcion regresa la serie de incrementos de una serie original
function increment_serie(serie::Array{Float64,1})
    inc = zeros(Float64,length(serie))
    for i = 1:(length(serie)-1)
        inc[i+1] = serie[i+1] - serie[i]
    end
    return inc
end
###############################################################################
#funcion que regresa la serie de signos de una serie original
function sign_serie(serie::Array{Float64,1})
    sgn = zeros(Float64,length(serie))
    for i = 1:(length(serie))
        sgn[i] = sign(serie[i])
    end
    return sgn - mean(sgn)
end
######################################################################################
#la funcion convierte la serie de entrada en su serie diferenciada
function differ_serie(serie::Array{Float64})
    n = length(serie)
    dif = zeros(Float64, n)
    dif[1] = serie[2] - serie[1] #se hacen los extremos
    dif[n] = (serie[n] - serie[n-1])
    for i = 2:(n-1)
        dif[i] = (serie[i+1] - serie[i-1]) / 2
    end
    return dif
end
###################################################################################
function polyfits(x,y,n) #funcion que ajusta un polinomio de grado n
  A = [ float(x[i])^p for i = 1:length(x), p = 0:n ]
  A \ y
end
#####################################################################################
function polyeval2(a,b,c,d,e,x) #funcion que evalua un polinomio de cuarto orden en x
  x^4 * a + x^3 * b + x^2 * c + x * d + e
end
#####################################################################################
function eval2(a,b,c,x)
    x^2 * a + x * b + c
end

########################################################################################################################################
#esta funcion hace un runnnig average de tamanio s en la serie, la ventana tiene
#que ser un numero impar
function run_avg(serie::Array{Float64}, s::Int64)
    n = length(serie)
    k = (s - 1) / 2
    x = n - 2 * k
    ra = Array(Float64, x)
    for i = (k+1):(n - k)
        ra[i-k] = mean(serie[(i-k):(i+k)])
    end
end
#################################################################################################
#la siguiente funcion hace un tipo binning sobre los datos con ventanas de tamanio m
#regresa un arreglo de menor tamanio al original pero distribuido homogeneamente
function binning_dfa(serie::Array{Float64}, m::Int64, nu::Float64)
    n = length(serie) #el log10(6) es necesario porque es el minimo
    bin = Int64[]
    #println(nu)
    #println(bin)
    j = 1
    for i = 0:(m-1)
        tam = nu / 2 + 0.01
        x = 0
        while (i * nu + nu / 2) > (serie[j]) && j != n
            if abs(serie[j] - (i * nu)) < tam
                #println(i,'\t', j,'\t',tam,'\t', serie[j], '\t', i * nu)
                x = j
                tam = abs(serie[j] - (i * nu))
            end
            j += 1
        end
        if x != 0
            push!(bin,x)
        end
    end
    #println(bin)
    return bin
end
#########################################################################
#the next function calculates the intervals in a time series Δx
function inter_serie(s::Array{Float64},d::Int64)
    n = length(s)-d
    dif = [99999.0 for i=1:n]
    for i = 1:n
        k = 1
        if s[i] != 0 && s[i+d] != 0
            dif[i] = s[i+1] - s[i]
        end
        if s[i] != 0 && s[i+d] == 0
            while s[i+k] == 0 && (i+k) < n
                k+=1
            end
            if s[i+k+d] != 0
                dif[i+k-1] = s[i+k+d] - s[i]
            end
        end
    end
    return dif
end
##########################################################################
#la siguiente funcion es solo para generar un arreglo que se pueda imprimir
#toma como entrada un diccionario.
function out_dic(d::Dict{Float64,Float64}, lb::Int64)
    if haskey(d, 99999.0)
        a = zeros(length(d)-1,3)
    else
        a = zeros(length(d),3)
    end
    k = 1
    for i in keys(d)
        if i == 99999.0; continue; end;
        a[k,1] = lb; a[k,2] = i; a[k,3] = d[i]
        k += 1
    end
    return a
end
