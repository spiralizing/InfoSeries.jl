#this is the main file
module InfoSeries
    using Clustering, DataFrames, Distances, StatsBase

    export cross_over,dif_ruidor,slopes,trasup_serie!,tras_serie!,integrate_serie!, integrate, magnitude_serie, increment_serie, sign_serie, differ_serie, polyfits, polyeval2, eval2, mv_avg, binning_dfa, inter_serie, out_dic, dfamod_calc, fast_fourier
    export rank_seriestxt,rank_freqtxt,rank_freq_cum,freq,rank_freq,prob_marg, prob_conj, prob_conj3, prob_claplace, aut_corr, cross_corr, mutual_inf, auto_inf, histo, distro, transfer, rms, find_line, find_crov
    export rec_swmap,rec_map,admat_out,init_swmat,get_slinks,get_glinks,vmin_blocks,vr_blocks, findgroup, v_blocks,visibility,inv_serie,hvi_blocks,freq_series,rank_series_i,rf_blocks,hv_links,hgroup_blocks,group_blocks,hv_blocks,reduce_hv,ham_distexp,isinarray,ham_dist,hd_visibility,sh_ent,red_hv,rank_series,inter_symbol,csvtointervs,h_visibility,R_squared,interval_serie, note_serie, rand_series,consonance_series,gaps_notas!, gaps_silencios!, min_voces, max_tempo, serie_notas, indice, notas_hertz!, series_desdoble, filt_vozcsv, nota_intens!, note_vec, dfa_calc, rounding!, csvtoserie
    export mat_trans1, mat_trans, symb_sh_ent, txt_series, tmat_txt, txt_words,markov_txt,random_spacing

    include("func_probs.jl")
    include("func_series.jl")
    include("constr_series.jl")
    include("func_stochastic.jl")
end
