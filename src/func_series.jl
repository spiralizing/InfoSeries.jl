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
################################################################################################################
#La siguiente funcion calcula el DFA, el primer argumento son las variables, el segundo es el eje del tiempo
#el tercer argumento es una opcion p=1 para la serie original, p=0 para las series de incremento y signo
function dfa_calc(notas::Array{Float64,2}, temp::Array{Float64,1}, p::Int64)
    tam = size(notas)
    serie = zeros(tam[1], tam[2])
    #se integra la serie
    #println(tam[2])
    #println("Passed! flag 0")
    for i = 1:tam[2]
        serie[:,i] = integrate(notas[:,i])
    end
    #println("Passed! flag 1")
    graficas = Array{Float64, 2}
    n = 4 #es el numero minimo de cajas
    tmax = convert(Int,trunc(tam[1] /n))
    tmin = 6
    si = zeros(tmax-tmin+1)
    fu = zeros(tmax-tmin+1)
    N = tam[1]
    for s = tmin:tmax
        nu = convert(Int,trunc(tam[1] / s)) - 1 #se calcula cuantos bloques de s elementos pueden obtenerse
        Fn = zeros(tam[1],tam[2])
        for x = 1:tam[2]
            for k = 0:nu
                block = serie[k * s + 1: k * s + s , x] #asigna los bloques de tamanio s al arreglo block
                blockt = temp[k * s + 1: k* s + s]
                coef1 = polyfits(blockt, block, 2) #se ajusta un polinomio al bloque
                for j = 1:s #se resta el polinomio al bloque
                    Fn[k * s + j , x] = (block[j] - eval2(coef1[3],coef1[2],coef1[1],blockt[j]))
                end
            end
        end
        Fs = 0
        for k = 1:N
            Vn = Fn' #se transpone la matriz Fn para convertirla en un vector y asi usar la funcion dot
            Fs += dot(Vn[1:end , k],Vn[1:end , k]) #se realiza el producto punto y se suma sobre todas las entradas
        end
        if p == 1
            Fnn = sqrt(Fs/ N) #se calcula la funcion F(s)
        else
            Fnn = sqrt(Fs / N) / N
        end
        si[s-tmin+1] = log10(s)
        fu[s-tmin+1] = log10(Fnn)
        #println(s ," ", Fnn)
    end
    #println("Paseed! flag 2")
    #writedlm("grafica.dat", [si fu], '\t')
    sid = si - log10(tmin)
    tm = round(Int,trunc(sid[end] / 0.05)) #0.08 es la ventana mas grande de la serie.
    ind = binning_dfa(sid, tm, 0.05)
    #println(length(ind))
    #println(ind[end])
    #println(length(fu))
    #println(length(si)"\t" length(fu))
    fy = zeros(length(ind))
    sx = zeros(length(ind))
    for z = 1:length(ind) #estos son los arreglos finales
        if ind[z] > ind[end] || ind[z] == 0 || ind[z] < 0; continue; end
        fy[z] = fu[ind[z]]
        sx[z] = si[ind[z]]
    end
    #println("Passed! flag 3")
    graficas = [sx fy] #la grafica que se exporta
    lin = polyfits(sx[1:end],fy[1:end],1) #ajuste lineal a los puntos
    return lin[2], graficas
end
##################################################################################################################################
#Para una sola voz
function dfa_calc(notas::Array{Float64,1}, temp::Array{Float64,1}, p::Int64)
    tam = length(notas)
    serie = zeros(tam)
    #se integra la serie
    #println(tam[2])
    serie = integrate(notas)
    graficas = Array{Float64, 2}
    n = 4 #es el numero minimo de cajas
    tmax = convert(Int,trunc(tam /n))
    tmin = 6
    si = zeros(tmax-tmin+1)
    fu = zeros(tmax-tmin+1)
    N = tam
    #println("Passed flag 0 dfa1")
    for s = tmin:tmax
        nu = convert(Int,trunc(tam[1] / s)) - 1 #se calcula cuantos bloques de s elementos pueden obtenerse
        Fn = zeros(tam)
        for k = 0:nu
            block = serie[k * s + 1: k * s + s] #asigna los bloques de tamanio s al arreglo block
            blockt = temp[k * s + 1: k* s + s]
            coef1 = polyfits(blockt, block, 2) #se ajusta un polinomio al bloque
            for j = 1:s #se resta el polinomio al bloque
                Fn[k * s + j] = (block[j] - eval2(coef1[3],coef1[2],coef1[1],blockt[j]))
            end
        end
        Fs = 0
        for k = 1:N
            Vn = Fn' #se transpone la matriz Fn para convertirla en un vector y asi usar la funcion dot
            Fs += dot(Vn[1:end , k],Vn[1:end , k]) #se realiza el producto punto y se suma sobre todas las entradas
        end
        if p == 1
            Fnn = sqrt(Fs/ N) #se calcula la funcion F(s)
        else
            Fnn = sqrt(Fs / N) / N
        end
        si[s-tmin+1] = log10(s)
        fu[s-tmin+1] = log10(Fnn)
        #println(s ," ", Fnn)
    end
    #writedlm("grafica.dat", [si fu], '\t')
    sid = si - log10(tmin)
    tm = round(Int,trunc(sid[end] / 0.05)) #0.08 es la ventana mas grande de la serie.
    ind = binning_dfa(sid, tm, 0.05)
    #println(length(ind))
    #println(ind[end])
    #println(length(fu))
    #println(length(si)"\t" length(fu))
    fy = zeros(length(ind))
    sx = zeros(length(ind))
    for z = 1:length(ind) #estos son los arreglos finales
        if ind[z] > ind[end] || ind[z] == 0 || ind[z] < 0; continue; end
        fy[z] = fu[ind[z]]
        sx[z] = si[ind[z]]
    end
    graficas = [sx fy] #la grafica que se exporta
    lin = polyfits(sx[1:end],fy[1:end],1) #ajuste lineal a los puntos
    return lin[2], graficas
end
#############################################################################################################
#El siguiente programa calcula el dfa modificado (signos y magnitudes)
function dfamod_calc(notas::Array{Float64,2}, temp::Array{Float64,1})
    tam = size(notas)
    incserie = zeros(tam[1],tam[2])
    sgnserie = zeros(tam[1],tam[2])
    magserie = zeros(tam[1],tam[2])
    graficas = Array(Array{Float64, 2}, 2) #es el arreglo de graficas
    #se necesita construir primero la serie descompuesta de magnitudes y signos
    #para esto, el metodo de ashkenazy sugiere primero asegurar que la serie
    #de incrementos este anticorrelacionada (α < 0.5)
    #se construye la serie de incrementos ΔX
    for k = 1:tam[2]
        incserie[:,k] = increment_serie(notas[:,k])
    end
    #se calcula el nDFA-2 de la serie de incrementos
    # si α > 0.5 se diferencia la serie hasta que α < 0.5
    al = dfa_calc(incserie, temp, 1)
    while al[1] > 0.5

        for k = 1:tam[2]
    	    incserie[:,k] = differ_serie(incserie[:,k])

        end

        al = dfa_calc(incserie, temp, 1)
    end
    #una vez tenida la serie de incrementos anticorrelacionada
    #se descompone en magnitudes y signos
    for k = 1:tam[2]
        magserie[:,k] = magnitude_serie(incserie[:,k])
        sgnserie[:,k] = sign_serie(incserie[:,k])
    end

    graficas[1] = dfa_calc(magserie, temp, 0)[2]
    graficas[2] = dfa_calc(sgnserie, temp, 0)[2]
    return graficas
end
#El siguiente programa calcula el dfa modificado (signos y magnitudes)
function dfamod_calc(notas::Array{Float64,1}, temp::Array{Float64,1})
    tam = length(notas)
    incserie = zeros(tam)
    sgnserie = zeros(tam)
    magserie = zeros(tam)
    graficas = Array(Array{Float64, 2}, 2) #es el arreglo de graficas
    #se necesita construir primero la serie descompuesta de magnitudes y signos
    #para esto, el metodo de ashkenazy sugiere primero asegurar que la serie
    #de incrementos este anticorrelacionada (α < 0.5)
    #se construye la serie de incrementos ΔX
    incserie = increment_serie(notas)
    #se calcula el nDFA-2 de la serie de incrementos
    # si α > 0.5 se diferencia la serie hasta que α < 0.5
    al = dfa_calc(incserie, temp, 1)
    while al[1] > 0.5
        incserie = differ_serie(incserie)
        al = dfa_calc(incserie, temp, 1)
    end
    #una vez tenida la serie de incrementos anticorrelacionada
    #se descompone en magnitudes y signos
    magserie = magnitude_serie(incserie)
    sgnserie = sign_serie(incserie)

    graficas[1] = dfa_calc(magserie, temp, 0)[2]
    graficas[2] = dfa_calc(sgnserie, temp, 0)[2]
    return graficas
end
################################################################################################################################
#root mean squared 1
function rms(f::Array{Float64,1}, p::Array{Float64,1})
    n = length(f)
    s = 0
    for i = 1:n
        s += (f[i] - p[i]) ^ 2 / n
    end
    return sqrt(s)
end

###################################################################################################################################
#La siguiente funcion calcula el numero de puntos en el cual hay un menor error en la desviacion estandar, d es el indice donde inicia
function find_line(puntos::Array{Float64,2}, d::Int64)
    n = size(puntos)[1]
    rs = ones(n)
    for i = d:n
        f = zeros(i)
        l = polyfits(puntos[1:i,1], puntos[1:i,2], 1)
        for j = 1:i
            f[j] = puntos[j,1] * l[2] + l[1]
        end
        rs[i] = rms(f, puntos[:,2])
    end
    return minimum(rs), indmin(rs)
end
######################################################################################################
#La siguiente funcion es para identificar los crossover que hay en el dfa, tam = tamanio - 1 de la ventana, numero par!
function find_crov(puntos::Array{Float64,2}, tam)
    x = []
    y = []
    m = []
    z = []
    n = size(puntos)[1]
    i = 1
    while (i + tam) < n
        f = zeros(tam+1)
        f2 = zeros(tam+3)
        l = polyfits(puntos[i:(i+tam),1], puntos[i:(i+tam),2], 1)
        for j = 1:(tam+1)
            f[j] = puntos[i+j-1] * l[2] + l[1]
        end
        push!(x, i + tam/2); push!(y, rms(f, puntos[i:(i+tam),2])); push!(m, l[2]);
        if i == 1
            for j = 1:(tam+3)
                f2[j] = puntos[i+j-1] * l[2] + l[1]
            end
            push!(z,rms(f2, puntos[i:(i+tam+2),2]))
        elseif (i+tam) < n-1
            for j = 1:(tam+3)
                f2[j] = puntos[i+j-2] * l[2] + l[1]
            end
            push!(z, rms(f2, puntos[(i-1):(i+tam+1),2]))
        else
            for j = 1:(tam+3)
                f2[j] = puntos[i+j-3] * l[2] + l[1]
            end
            push!(z,rms(f2, puntos[(i-2):(i+tam),2]))
        end
        #println(i + tam/2,'\t', rms(f, puntos[i:(i+tam),2]))
        i += 1
    end
    mn = zeros(length(m))
    for i = 1:length(m)-1
        mn[i+1] = abs(m[i+1] - m[i])/ abs(m[i])
    end
    return [x y z mn]
end
###########################################################################################################
function fast_fourier(s::Array{Float64,1})
    ind = log10(1:length(s))
    tam = length(s)
    f = log10(abs2(fft(s)))
    n = Int64(ceil(length(f) / 2))
    tm = floor(Int64, (ind[n]-ind[10])/ 0.03)
    idx = append!(collect(1:11),binning_dfa(ind[11:n], tm, 0.05) + 11)
    fo = Array(Float64,length(idx))
    indo = Array(Float64,length(idx))
    for i = 1:length(idx)
        fo[i] = f[idx[i]]
        indo[i] = ind[idx[i]]
    end
    return [indo fo]
end
#################################################################################################################
function tras_serie!(s::Array{Float64,2})
    mi = zeros(size(s)[2])
    for i = 1:(size(s)[2]-1)
        mi[i] = minimum(filter(x -> x!=0, s[:,i+1]))-1
    end
    for i = 1:size(s)[1]
        for j = 2:size(s)[2]
            if s[i,j] == 0; continue; end
            s[i,j] = s[i,j]-mi[j-1]
        end
    end
    return s
end
########################################################################################
function trasup_serie!(s::Array{Float64,2},s2::Array{Float64,2}) #mueve una serie a la altura original, se toma como referencia la serie original
    mi = zeros(size(s)[2])
    for i = 1:(size(s)[2]-1)
        mi[i] = minimum(filter(x -> x!=0, s2[:,i+1]))-1
    end
    for i = 1:size(s)[1]
        for j = 2:size(s)[2]
            if s[i,j] == 0; continue; end
            s[i,j] = s[i,j]+mi[j-1]
        end
    end
    return s
end
######################################################################################
function slopes(s::Array{Float64,2}) # calcula la serie de pendientes a puntos vecinos.
    n = size(s)[1]
    m = zeros(n)
    m[1] = (s[2,2] - s[1,2])/(s[2,1] - s[1,1]) #se hacen los extremos
    m[n] = (s[n,2] - s[n-1,2])/(s[n,1] - s[n-1,1])
    for i = 2:(n-1)
        m[i] = (s[i+1,2] - s[i,2])/(s[i+1,1] - s[i,1])
    end
    return m
end
#####################################################################################################################
#La siguiente funcion hace una simple resta entre la funcion de fluctuaciones del DFA F(s) y la recta con pendiente 0.5 correspondiente a ruido descorrelacionado.
function dif_ruidor(s::Array{Float64,2})
    n = size(s)[1]
    dif = zeros(n)
    b = s[1,2] - s[1,1] * 0.5
    for i = 1:n
        dif[i] = s[i,2] - (0.5 * s[i,1] + b)
    end
    return dif
end
####################################################################################################
function rand_series(s::Array{Float64,2}) #este programa hace un random shuffle sobre cada una de las voces de una serie.
    out = zeros(size(s)[1],size(s)[2])
    out[:,1] = s[:,1]
    for i in 2:size(s)[2]
        out[:,i] = shuffle(s[:,i])
    end
    return out
end
