#this is the main file
module InfoSeries
    
    export integrate_serie!, integrate, magnitude_serie, increment_serie, sign_serie, differ_serie, polyfits, polyeval2, eval2, run_avg, binning_dfa, inter_serie, out_dic
    export prob_marg, prob_conj, prob_claplace, aut_corr, cross_corr, mutual_inf, auto_inf, histo

    include("func_probs.jl")
    include("func_series.jl")
end
