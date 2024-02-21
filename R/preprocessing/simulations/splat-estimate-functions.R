library(tidyverse)
library(magrittr)
library(splatterBatch)
library(Matrix)
library(edgeR)

winsorize = function(x, q) {
	checkmate::check_numeric(x, any.missing = FALSE)
	checkmate::check_number(q, lower = 0, upper = 1)
	lohi = stats::quantile(x, c(q, 1 - q), na.rm = TRUE)
	if (diff(lohi) < 0) { lohi = rev(lohi) }
	x[!is.na(x) & x < lohi[1]] = lohi[1]
	x[!is.na(x) & x > lohi[2]] = lohi[2]
	return(x)
}

logistic = function(x, x0, k) {
	1 / (1 + exp(-k * (x - x0)))
}

splatEstMean = function(norm.counts, params) {

	means = rowMeans(norm.counts)
	means = means[means != 0]

	means = winsorize(means, q = 0.1)

	fit = fitdistrplus::fitdist(means, "gamma", method = "mge",
								 gof = "CvM")
	if (fit$convergence > 0) {
		warning("Fitting means using the Goodness of Fit method failed, ",
				"using the Method of Moments instead")
		fit = fitdistrplus::fitdist(means, "gamma", method = "mme")
	}

	params = setParams(params, mean.shape = unname(fit$estimate["shape"]),
						mean.rate = unname(fit$estimate["rate"]))

	return(params)
}

splatEstLib = function(counts, params) {

	lib.sizes = colSums(counts)

	if (length(lib.sizes) > 5000) {
		message("NOTE: More than 5000 cells provided. ",
				"5000 sampled library sizes will be used to test normality.")
		lib.sizes.sampled = sample(lib.sizes, 5000, replace = FALSE)
	} else {
		lib.sizes.sampled = lib.sizes
	}

	norm.test = shapiro.test(lib.sizes.sampled)
	lib.norm = norm.test$p.value > 0.2

	if (lib.norm) {
		fit = fitdistrplus::fitdist(lib.sizes, "norm")
		lib.loc = unname(fit$estimate["mean"])
		lib.scale = unname(fit$estimate["sd"])
		message("NOTE: Library sizes have been found to be normally ",
				"distributed instead of log-normal. You may want to check ",
				"this is correct.")
	} else {
		fit = fitdistrplus::fitdist(lib.sizes, "lnorm")
		lib.loc = unname(fit$estimate["meanlog"])
		lib.scale = unname(fit$estimate["sdlog"])
	}

	params = setParams(params, lib.loc = lib.loc, lib.scale = lib.scale,
						lib.norm = lib.norm)

	return(params)
}

splatEstOutlier = function(norm.counts, params) {

	means = rowMeans(norm.counts)
	lmeans = log(means)

	med = median(lmeans)
	mad = mad(lmeans)

	bound = med + 2 * mad

	outs = which(lmeans > bound)

	prob = length(outs) / nrow(norm.counts)

	params = setParams(params, out.prob = prob)

	if (length(outs) > 1) {
		facs = means[outs] / median(means)
		fit = fitdistrplus::fitdist(facs, "lnorm")

		params = setParams(params,
							out.facLoc = unname(fit$estimate["meanlog"]),
							out.facScale = unname(fit$estimate["sdlog"]))
	}

	return(params)
}

splatEstBCV = function(counts, params) {

	# Add dummy design matrix to avoid print statement
	design = matrix(1, ncol(counts_filtered), 1)
	disps = edgeR::estimateDisp(counts_filtered, design = design)

	params = setParams(params,
						bcv.common = 0.1 + 0.25 * disps$common.dispersion,
						bcv.df = disps$prior.df)

	return(params)
}

splatEstDropout = function(norm.counts, params) {

	means = rowMeans(norm.counts)

	x = log(means)

	obs.zeros = rowSums(norm.counts == 0)

	y = obs.zeros / ncol(norm.counts)

	df = data.frame(x, y)

	fit = nls(y ~ logistic(x, x0 = x0, k = k), data = df,
			   start = list(x0 = 0, k = -1))

	#exp.zeros = dnbinom(0, mu = means, size = 1 / 0.1) * ncol(norm.counts)

	#present = max(obs.zeros - exp.zeros) > 0.1 * ncol(norm.counts)

	mid = summary(fit)$coefficients["x0", "Estimate"]
	shape = summary(fit)$coefficients["k", "Estimate"]

	params = setParams(params, dropout.mid = mid, dropout.shape = shape)

	return(params)
}

splatEstimate.matrix = function(counts, params = newSplatParams()) {

	checkmate::assertClass(params, "SplatParams")

	# Normalise for library size and remove all zero genes
	lib.sizes = colSums(counts)
	lib.med = median(lib.sizes)
	norm.counts = t(t(counts) / lib.sizes * lib.med)
	norm.counts = norm.counts[rowSums(norm.counts > 0) > 1, ]

	params = splatEstMean(norm.counts, params)
	params = splatEstLib(counts, params)
	params = splatEstOutlier(norm.counts, params)
	params = splatEstBCV(counts, params)
	params = splatEstDropout(norm.counts, params)

	params = setParams(params, nGenes = nrow(counts),
						batchCells = ncol(counts))
	return(params)
}
