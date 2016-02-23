#aqui hay varias funciones de probabilidades y calculos de teoria de informacion
################################################################################
#Esta funcion calcula la probabilidad marginal de una variable aleatoria
#toma como entrada una serie y regresa un arreglo de marginales
function prob_marg(serie::Array{Float64,1})
    tam = length(serie)
    M = Dict{Float64,Float64}()
    for i = 1:tam
        M[serie[i]] = get(M, serie[i], 0) + 1
    end
    for j in keys(M)
        M[j] = M[j] / tam
    end
    return M
end
###############################################################################
#esta funcion calcula la probabilidad conjunta de dos variables en dos series
#, separadas por una distancia d > 0
function prob_conj(s1::Array{Float64}, s2::Array{Float64}, d::Int64)
    tam = min(length(s1), length(s2))
    Pxy = Dict{String, Float64}()
    Pyx = Dict{String, Float64}()
    for i = 1:tam-d
        xy = string(s1[i], ',', s2[i+d]) #se calculan las probabilidades conjuntas
        yx = string(s2[i], ',', s1[i+d])
        Pxy[xy] = get(Pxy, xy, 0) + 1/(tam-d)
        Pyx[yx] = get(Pyx, yx, 0) + 1/(tam-d)
    end
    return Pxy, Pyx
end
############################################################################
#esta funcion calcula la probabilidad conjunta con la regla de sucesor de laplace
# p(x,y) = (n + 1 )/ (N + c^2)
function prob_claplace(s1::Array{Float64}, s2::Array{Float64}, d::Int64)
    c = (length(prob_marg(s1)), length(prob_marg(s2)))
    tam = min(length(s1), length(s2))
    nxy = Dict{String, Float64}()
    nyx = Dict{String, Float64}()
    Pxy = Dict{String, Float64}()
    Pyx = Dict{String, Float64}()
    for i = 1:tam-d
        xy = string(s1[i], ',', s2[i+d]) #se calculan las probabilidades conjuntas
        yx = string(s2[i], ',', s1[i+d])
        nxy[xy] = get(nxy, xy, 0) + 1
        nyx[yx] = get(nyx, yx, 0) + 1
    end
    for j in keys(nxy)
        Pxy[j] = ( nxy[j] + 1 ) / ((tam - d) + (c[1] * c[2]))
    end
    for j in keys(nyx)
        Pyx[j] = ( nyx[j] + 1 ) / ((tam - d) + (c[1] * c[2]))
    end
    return Pxy, Pyx
end
#################################################################################
#this calculates the auto-correlation function in terms of a delay τ
function aut_corr(serie::Array{Float64})
    μ = mean(serie)
    σ = 0 #this is the variance (sigma squared)
    n = length(serie)
    ac = zeros(n,2)
    for i = 1:n
        σ += 1 / n * (serie[i] - μ) ^ 2
    end
    #println(σ)
    for d = 0:(n-1)
        suma = 0
        for i = 1:(n-d)
            suma += (serie[i] - μ) * (serie[i+d] - μ)
        end
        ac[d+1,1] = d
        ac[d+1,2] =  1 / (n-d) * suma / σ
    end
    return ac
end
################################################################################
#this is the cross correlation function
function cross_corr(s1::Array{Float64}, s2::Array{Float64})
    μ1 = mean(s1); μ2 = mean(s2)
    σ1 = 0; σ2 = 0
    n = length(s1)
    cc = zeros(2 * n - 1, 2)
    for i = 1:n
        σ1 += 1 / n * (s1[i] - μ1) ^ 2
        σ2 += 1 / n * (s2[i] - μ2) ^ 2
    end
    σ1 = σ1 ^ (1 / 2)
    σ2 = σ2 ^ (1 / 2)
    for d = 0:(n-1)
        sum1 = 0; sum2 = 0
        for i = 1:(n-d)
            sum1 += (s1[i] - μ1) * (s2[i+d] - μ2)
            sum2 += (s1[i+d] - μ1) * (s2[i] - μ2)
        end
        cc[n-d,1] = - d; cc[n-d,2] = 1 / (n-d) * sum2 / (σ1 * σ2)
        cc[n+d,1] = d; cc[n+d,2] = 1 / (n-d) * sum1 / (σ1 * σ2)
    end
    cc[n+1] = cc[n+1] / 2
    return cc
end
##################################################################################
#the next script calculates the mutual information function of a delay between 2 series
function mutual_inf(s1::Array{Float64}, s2::Array{Float64})
    tam = min(length(s1), length(s2))
    My = prob_marg(s2)
    Mx = prob_marg(s1)
    ind = zeros(2 * tam - 1)
    IM = zeros(2 * tam -1)
    for d = 0:(tam-1)
        Pxy = Dict{AbstractString,Float64}() #se definen los diccionarios para las probabilidades conjuntas
        Pyx = Dict{AbstractString,Float64}()
        for i = 1:(tam-d)
            xy = string(s1[i], ',', s2[i+d]) #se calculan las probabilidades conjuntas
            yx = string(s2[i], ',', s1[i+d])
            Pxy[xy] = get(Pxy, xy, 0) + 1/(tam-d)
            Pyx[yx] = get(Pyx, yx, 0) + 1/(tam-d)
        end
        for j in keys(Pxy) #se calcula la informacion mutua pa delante
            A = split(j, ',')
            x = parse(Float64,A[1])
            y = parse(Float64,A[2])
            IM[tam + d] += Pxy[j] * log10(Pxy[j]/(Mx[x] * My[y]))
        end
        for j in keys(Pyx) #se calcula la informacion mutua pa tras
            A = split(j, ',')
            x = parse(Float64,A[2])
            y = parse(Float64,A[1])
            IM[tam - d] += Pyx[j] * log10(Pyx[j]/(Mx[x] * My[y]))
        end
        ind[tam - d] = - d
        ind[tam + d] = d
    end
    IM[tam] = IM[tam] / 2
    return [ind IM]
end
##################################################################################
#This is the calculation of the auto mutual information.
function auto_inf(s::Array{Float64})
    tam = length(s)
    M = prob_marg(s)
    ind = zeros(tam)
    AI = zeros(tam)
    for d = 0:(tam-1)
        Pxx = Dict{AbstractString,Float64}() #se definen los diccionarios para las probabilidades conjuntas
        for i = 1:(tam-d)
            xx = string(s[i], ',', s[i+d]) #se calculan las probabilidades conjuntas
            Pxx[xx] = get(Pxx, xx, 0) + 1
        end
        for j in keys(Pxx)
            Pxx[j] = Pxx[j] / (tam - d)
        end
        for j in keys(Pxx) #se calcula la informacion mutua pa delante
            A = split(j, ',')
            x = parse(Float64,A[1])
            y = parse(Float64,A[2])
            AI[d+1] += Pxx[j] * log10(Pxx[j]/(M[x] * M[y]))
        end
        ind[d+1] = d
    end
    return [ind AI]
end
#################################################################################
#the next function is for the histogram of values for the series s
function histo(s::Array{Float64})
    H = Dict{Float64,Float64}()
    for i = 1:length(s)
        H[s[i]] = get(H, s[i], 0) + 1
    end
    return H
end
