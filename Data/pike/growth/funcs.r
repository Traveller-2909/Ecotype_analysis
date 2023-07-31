library(data.table)
library(gridExtra)
library(ggplot2)

mypallette <- viridis::viridis(6, direction = 1)

pl_theme =	theme_classic() +
  theme(panel.grid.major = element_line(),
        legend.box.background = element_rect(),
        axis.line = element_line(size = .5, colour = "black", linetype = 1),
        axis.ticks = element_line(size = .5, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5, face = "bold"),
        legend.key.size = unit(.5,"line"),
        legend.position = c(0.3,0.8))

init_fun = function() {
	list(
		LInf = runif(1, 2, 4),
		tZero = runif(1, 0.1, 1.5),
		k = runif(1, 0.05, 0.15),
		LInf_eco = runif(4, 2, 4),
		tZero_eco = runif(4, 0.1, 1.5),
		k_eco = runif(4, 0.05, 0.15),
		LInf_i = runif(120, 2, 4),
		tZero_i = runif(120, 0.1, 1.5),
		k_i = runif(120, 0.05, 0.15),
		LInf_phi = runif(1, 0, 2),
		tZero_phi =  runif(1, 0, 0.5),
		k_phi = runif(1, 0, 0.1),
		LInf_eco_phi = runif(4, 0.2, 2),
		tZero_eco_phi = runif(4, 0.2, 2),
		k_eco_phi = runif(4, 0, 0.1),
		phi = runif(1, 0, 0.1)
	)
}


Incpredict = function(pars, Agei) pars[1] * (1 - exp(-pars[2] * (Agei + pars[3])))

# x is a stanfit
ecotype_predict = function(x, quants = c(0.5, 0.05, 0.95), names = c("median", "lower", "upper")) {
	pars_eco = c("LInf_eco", "k_eco", "tZero_eco")
	samps = as.matrix(x, pars = pars_eco)

	# list of 6 ecotypes
	# each element is a matrix, 20000 rows (samples) by 15 age increments
	eco_pr = lapply(0:5, \(i) t(apply(samps[, c(1, 7, 13) + i], 1, Incpredict, Agei = 1:15)))

	# convert to a data table summarizing by quantiles
	eco_q = lapply(eco_pr, \(x) {
		res = data.frame(t(apply(x, 2, quantile, quants)))
		colnames(res) = names
		res$Agei = 1:15
		res
	})
	rbindlist(eco_q, idcol = "eco_id")
}

fish_predict = function(x, fish_lookup) {
	pars_fish = c("LInf_i", "k_i", "tZero_i")
	samps_f = as.matrix(x, pars = pars_fish)
	fish_pr = lapply(0:85, \(i) t(apply(samps_f[, c(1, 87, 173) + i], 1, Incpredict, Agei = 1:15)))
	fish_m = lapply(fish_pr, \(x) data.frame(fit = apply(x, 2, median), Agei = 1:15))
	fish_plot = rbindlist(fish_m, idcol = "fish_id")	
	fish_plot$eco_id = fish_lookup$ecotype[match(fish_plot$fish_id, fish_lookup$id3)]
	fish_plot
}

ecotype_curves = function(x, dat) {
	fish_lookup = unique(dat[, c("id3", "ecotype")])
	fish_lookup = fish_lookup[order(fish_lookup$id3),]
	eco_plot = ecotype_predict(x)
	fish_plot = fish_predict(x, fish_lookup)
	
	pl_curves = list()
	for(i in 1:6) {
		pl_curves[[i]] = ggplot(data = eco_plot[eco_id == i]) + 
				geom_ribbon(aes(x = Agei, ymin = lower, ymax = upper), fill = mypallette[i], alpha = 0.3) + 
				geom_line(aes(x = Agei, y = median), color = mypallette[i], linewidth = 2)
		pl_curves[[i]] = pl_curves[[i]] + 
			geom_line(data = fish_plot[eco_id == i], aes(x = Agei, y = fit, group = fish_id), color = mypallette[i], linewidth = 0.5) + 
			geom_point(data = dat[dat$ecotype == i, ], 
					   aes(x = agei,y = radmm, group = id), size = 1, alpha = .9, color="black") + 
			pl_theme + labs(x="Age", y="otolith increment (mm)")
	}
	pl_curves$nrow = 2
	do.call(grid.arrange, pl_curves)
}


