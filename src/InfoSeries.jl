#this is the main file
module InfoSeries

    export dif_ruidor,slopes,trasup_serie!,tras_serie!,integrate_serie!, integrate, magnitude_serie, increment_serie, sign_serie, differ_serie, polyfits, polyeval2, eval2, run_avg, binning_dfa, inter_serie, out_dic, dfamod_calc, fast_fourier
    export prob_marg, prob_conj, prob_conj3, prob_claplace, aut_corr, cross_corr, mutual_inf, auto_inf, histo, distro, transfer, rms, find_line, find_crov
    export gaps_notas!, gaps_silencios!, min_voces, max_tempo, serie_notas, indice, notas_hertz!, series_desdoble, filt_vozcsv, nota_intens!, note_vec, dfa_calc, rounding!, csvtoserie
    export mat_trans1, mat_trans

    include("func_probs.jl")
    include("func_series.jl")
    include("constr_series.jl")
    include("func_stochastic.jl")
end
