library(diffr)
filename1 = "/Users/jiancao/Desktop/new2/TFDTW.R"
filename2 = "/Users/jiancao/Desktop/old/TFDTW.R"
diffr(filename1, filename2)

filename1 = "/Users/jiancao/Desktop/new/compare.R"
filename2 = "/Users/jiancao/Desktop/old/compare.R"
diffr(filename1, filename2)


unit = as.character(item$unit[1])
x = item$value
x_raw = item$value_raw
y = y_processed

res = TwoStepDTW(x, y, t_treat, k, n_dtw1,
                 normalize_method = normalize_method,
                 ma = ma, ma_na = ma_na,
                 type_dtw1 = type_dtw1,
                 dist_quantile = dist_quantile,
                 step.pattern1 = step.pattern1, 
                 step.pattern2 = step.pattern2, 
                 plot_figures = plot_figures, ...)

x_warped = c(warp_using_weight(x_raw[1:res$cutoff], res$weight_a)[1:t_treat],
             warp_using_weight(x_raw[-(1:(res$cutoff - 1))],
                               res$avg_weight)[-1])


print(paste0("x="))
print(x)
print(paste0("y="))
print(y)
print(paste0("k=", k))
print(paste0("n_dtw1=", n_dtw1))
print(paste0("t_treat=", t_treat))
print(paste0("normalize_method=", normalize_method))
print(step.pattern)
print(paste0("plot_figures=", plot_figures))


