#' @title Extract Model Coefficients From a Fitted Mixture of Experts Model
#'
#' @description Extract model coefficients from a fitted mixture of experts model.
#'
#' @param object [\code{\link{EMglmnet}}]\cr
#'  An object of class \code{EMglmnet} of which to extract the model coefficients.
#' @param ... [\code{ANY}]\cr
#'  Arguments to \code{\link[glmnet]{coef.glmnet}}.
#'
#' @return A list with model coefficients belonging to the gating and expert models.
#'
#' @export

coef.EMglmnet = function(object, ...) {
	cfs = vector(mode = "list", length = 2)
	names(cfs) = c("gating", "experts")
	if (!is.null(object$gating))
		cfs$gating = coef(object$gating, ...)	
	cfs$experts = lapply(object$experts, coef, ...)
	cfs
}



# @title Plot Coefficients of a Fitted Mixture of Experts Model
#
# @description Plot coefficients of a fitted mixture of experts model.
#
# @param x [\code{\link{EMglmnet}}]\cr
#  An object of class \code{EMglmnet} of which to extract the model coefficients.
# @param ... [\code{ANY}]\cr
#  Arguments to \code{\link[glmnet]{coef.EMglmnet}}.
#
# @export

# plotCoefs = function(x, ...) {
	# cfs = coef(x, ...)
# ## lattice plot mit fixen coefs?
	
# }



#' @title Path Plots for a Fitted Mixture of Experts Model
#'
#' @description Path plots for a fitted mixture of experts model.
#'
#' @param x [\code{\link{EMglmnet}}]\cr
#'  An object of class \code{EMglmnet} of which to extract the model coefficients.
#' @param ... [\code{ANY}]\cr
#'  Arguments to \code{\link[glmnet]{plot.glmnet}}.
#'
#' @export

plotPath = function(x, ...) {
	## gating
	if (!is.null(x$gatingPath)) {
		if (x$J == 2)
			plot(x$gatingPath, ...)
		else {
			par(mfrow = c(2, ceiling(x$J/2)))
			plot(x$gatingPath, ...)
		}
	}
	## experts
	sapply(x$expertsPath, function(z) {
		if (x$K == 2)
			plot(z, ...)
		else {
			par(mfrow = c(2, ceiling(x$K/2)))
			plot(z, ...)
		}
	})
}


## plotPathInteractions


## extract single active variables?

#' @title Extract Interactions From a Fitted Mixture of Experts Model
#'
#' @description Extract interactions from a fitted mixture of experts model.
#'
#' @param object [\code{\link{EMglmnet}}]\cr
#'  An object of class \code{EMglmnet} of which to extract the model coefficients.
#' @param ... [\code{ANY}]\cr
#'  Arguments to \code{\link[glmnet]{coef.glmnet}}.
#'
#' @return A \code{data.frame} of active interactions.
#  list with model coefficients belonging to the gating and expert models.
#'
#' @export

getInteractions = function(object, ...) {
	cfs = coef(object, ...)
## FIXME: auf Standardisierung Rücksicht nehmen?
## FIXME: wegen Identifizierbarkeit transformieren?
	if (object$J == 2) {
		if (object$K == 2) {
			int = getInteractionsHelper(1:2)
		} else if (object$K > 2) {
			stop("not yet")
		} else {
			stop("something is wrong with 'object$K'")
		}
	} else if (object$J > 2) {
		if (object$K == 2) {
			combs = combn(object$J, 2)
			int = apply(combs, 2, getInteractionsHelper)
			# do.call rbind, rbind.fill
		} else if (object$K > 2) {
			stop("not yet")
		} else {
			stop("something is wrong with 'object$K'")
		}
	} else {
		stop("something is wrong with 'object$J'")
	}
	int	
}



#' @noRd

## FIXME: odds ratios
## was bedeutet softmax?
## j-te gating trennt was?
getInteractionsHelper = function(z) {
	diffExperts = cfs$experts[[z[2]]] - cfs$experts[[z[1]]]
	diffGating = cfs$gating[[z[2]]] - cfs$gating[[z[1]]]
	indExperts = glmnet::nonzeroCoef(diffExperts)
	indExperts = indExperts[indExperts != 1]
	indGating = glmnet::nonzeroCoef(diffGating)
	indGating = indGating[indGating != 1]
# print(indGating)
# print(indExperts)
	if (!length(indGating) | !length(indExperts))
		int = data.frame()
	else {
		int = cbind(
			expand.grid(gatingVar = rownames(diffGating)[indGating], expertVar = rownames(diffExperts)[indExperts]),
			expand.grid(gating = diffGating[indGating], experts = diffExperts[indExperts])
		)
		ord = order(abs(int$gating * int$experts), decreasing = TRUE)
	}
	int
}


### Pareto plot 

# plotInteractions(x, ...) {
	# int = getInteractions(x)
	# plot(exp(gating) ~ exp(experts), int)
	# ## labels, Betrag?
	
# }


## joint odds ratio of two variables



	
