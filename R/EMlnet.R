## FIXME: parameter that gives ratio between gating and experts lambdas
## FIXME: multinomial: CEM steps? Ng & McLachlan, Peralta & Soto


#' @title Penalized Mixtures of Experts
#'
#' @description
#' Fits a penalized mixture of experts model using the EM algorithm.
#' In the EM steps penalized linear models are fitted using function \code{\link[glmnet]{glmnet}} from \pkg{glmnet} is used.
#'
#' @param X A \code{matrix} or \code{data.frame} or \code{Matrix} containing
#'   the explanatory variables (must not contain an intercept column).
#' @param y A \code{factor} specifying the class membership for each observation.
#' @param colsGating Names or indices of columns in \code{X} to be used for the gating model. Default is all columns.
#' @param colsExperts Names or indices of columns in \code{X} to be used for the expert models. Default is all columns.
#' @param interceptGating Logical. Does the gating model include an intercept?
#'   Defaults to \code{TRUE}.
#' @param interceptExperts Logical. Does the expert model include an intercept?
#'   Defaults to \code{TRUE}.
#' @param offsetGating Offset term for the gating model.
#   A numeric vector. FIXME: matrix?
#' @param offsetExperts Offset term for the expert model.
#   A numeric vector. FIXME: matrix?
#' @param J The number of experts / mixture components. Defaults to 2.
#' @param lambda Penalty parameter. Can be a scalar or a vector of length \code{J+1} with different components for
#'   the \code{J} experts and the gating model. All components must be >= 0.
#' @param alpha Mixing parameter for the elastic net penalty. Can be a scalar or a vector of length \code{J+1} with different
#'   components for the \code{J} experts and the gating model. All components must be in [0,1].
#'   Defaults to 1.
#' @param standardize Logical. Should the columns of \code{X} be standardized prior to fitting the model? 
#   when necessary? return on original scale? macht das sinn fuer spase X?
#' @param iter.max Maximum number of EM iterations.
#' @param tol Small positive value to detect convergence of the EM algorithm.
#' @param init.weights An \code{nrow(X)} times \code{J} matrix containing initial group membership values for the EM algorithm.
#'  The row sums must be 1. Defaults to \code{NULL}, which means that initial weights will be obtained by random segmentation.
#'
#' @return An object of class \code{EMlnet}.
#'  A \code{list} with entries:
#'   \item{gating}{The gating model. An object of class \code{glmnet}.}
#'   \item{experts}{A \code{list} of \code{J} expert models. These are objects of class \code{glmnet}}
#'   \item{weights}{The final group membership values.}
#'   \item{iters}{The number of iterations in the EM algorithm.}
#'   \item{J,K, colsGating, colsExperts, lambda, offsetGating, offsetExperts}{See arguments.}
#'   \item{lev1}{Class labels present in the data.}
#'   \item{lev}{Class labels.}
#'
#' @export

EMlnet = function(y, X, colsGating = 1:ncol(X), colsExperts = 1:ncol(X), J = 2, lambda, alpha = 1, standardize = FALSE, interceptGating = TRUE, interceptExperts = TRUE, offsetGating = NULL, offsetExperts = NULL, iter.max = 100, tol = 1e-06, init.weights = NULL) {
	
	## check X
    if (is.null(dim(X))) 
        stop("'X' is not a matrix")
	if (!is.matrix(X) & !inherits(X, "Matrix"))
		X = as.matrix(X)
	N = nrow(X)

	## check y
    if (N != length(y))
        stop("'nrow(X)' and 'length(y)' are different")
    if (!is.factor(y)) 
        warning("'y' was coerced to a factor")
	y = as.factor(y)
    lev = lev1 = levels(y)
    counts = as.vector(table(y))
    if (any(counts == 0)) {
        empty = lev[counts == 0]
        warning(sprintf(ngettext(length(empty), "group %s is empty", 
            "groups %s are empty"), paste(empty, collapse = ", ")), 
            domain = NA)
        lev1 = lev[counts > 0]
        y = factor(y, levels = lev1)
        counts = as.vector(table(y))
    }
    if (length(lev1) < 2L) 
        stop("training data from only one group given")
    names(counts) = lev1
    prior = counts/N
	K = length(lev1)
	## check for groups with small numbers of training observations
	if (any(counts == 1))
		stop("One or more class has only 1 training observation")
	
	## class indicator vectors
	Y = diag(K)[y,]
	colnames(Y) = lev1
	Y[Y == 0] = .Machine$double.eps
	Y = log(Y)
	## in the binary case keep only the higher factor level	
	if (K == 2) {
		Y = Y[,2]
	}

	## lambda
	if (any(lambda < 0))
		stop("'lambda' must be > 0")
	if (length(lambda) == 1) {
		lambda = rep(lambda, 1 + J)
	} else if (length(lambda) != 1 + J) {
		stop("'length(lambda)' is not 1 + J")
	}

	## alpha
	if (any(alpha < 0 | alpha > 1))
		stop("'alpha' must be in [0,1]")
	if (length(alpha) == 1) {
		alpha = rep(alpha, 1 + J)
	} else if (length(alpha) != 1 + J) {
		stop("'length(alpha)' is not 1 + J")
	}

	## initialization
	if (!is.null(init.weights)) {
		## check weights
		if (nrow(init.weights) != N)
			stop("'nrow(init.weights)' and 'nrow(X)' must be equal")
		if (ncol(init.weights) != J)
			stop("'ncol(init.weights)' must equal 'J'")
		s = rowSums(init.weights)
		if (any(s != 1)) {# all.equal?, tolerance?
			warning("'init.weights' is rescaled such that all row sums are 1")
			init.weights = init.weights/s
		}
		weights = init.weights
	} else {
		## random segmentation
		weights = matrix(runif(N*J), ncol = J)
		weights = weights/rowSums(weights)
	}
	
	weight.targets = weights
	weight.targets[weights < .Machine$double.eps] = .Machine$double.eps ## literature says 1e-03 works well
	## in the binary case keep only the second group
	weight.targets = log(weight.targets)
	if (J == 2) {
		weight.targets = weight.targets[,2]
	}

	oldPars = newPars = 0

	for (i in 1:iter.max) {
		## M step
		## calculate models for lambda sequence, extract coefficients and fitted values for the user-defined lambda
		if (J == 2) {
			gating = glmnet::glmnet(X[,colsGating, drop = FALSE], weight.targets, family = "gaussian", alpha = alpha[J+1],
				intercept = interceptGating, offset = offsetGating, standardize = standardize, lambda = lambda[J+1])
## FIXME: max und min lambda
			if (interceptGating)
				parsGating = c(gating$a0, as.vector(gating$beta))
			else
				parsGating = as.vector(gating$beta)
		} else {
			gating = parsGating = list()
			for (j in 1:J) {
				gating[[j]] = glmnet::glmnet(X[,colsGating, drop = FALSE], weight.targets[,j], family = "gaussian",
					alpha = alpha[J+1], intercept = interceptGating, offset = offsetGating, standardize = standardize,
					lambda = lambda[J+1])
				if (interceptGating)
					parsGating[[j]] = c(as.numeric(gating[[j]]$a0), as.numeric(sapply(gating[[j]]$beta, as.vector)))
				else
					parsGating[[j]] = as.numeric(sapply(gating$beta, as.vector))
			}
			parsGating = unlist(parsGating)
		}
		experts = parsExperts = list()
		if (K == 2) {
			for (j in 1:J) {
				experts[[j]] = glmnet::glmnet(X[,colsExperts, drop = FALSE], Y, weights = weights[,j], family = "gaussian",
					alpha = alpha[j], intercept = interceptExperts, offset = offsetExperts, standardize = standardize,
					lambda = lambda[j])
				if (interceptExperts)
					parsExperts[[j]] = c(experts[[j]]$a0, as.vector(experts[[j]]$beta))
				else 
					parsExperts[[j]] = as.vector(experts[[j]]$beta)
			}
		} else {
			for (k in 1:K) {
				experts[[k]] = parsExperts[[k]] = list()
				for (j in 1:J) {
					experts[[k]][[j]] = glmnet::glmnet(X[,colsExperts, drop = FALSE], Y[,k], weights = weights[,j],
						family = "gaussian", alpha = alpha[j], intercept = interceptExperts, offset = offsetExperts,
						standardize = standardize, lambda = lambda[j])
					if (interceptExperts)
						parsExperts[[k]][[j]] = c(as.numeric(experts[[k]][[j]]$a0), as.numeric(sapply(experts[[k]][[j]]$beta, as.vector)))
					else 
						parsExperts[[k]][[j]] = as.numeric(sapply(experts[[k]][[j]]$beta, as.vector))
				}
			}
		}
		parsExperts = unlist(parsExperts)

		## termination criterion
		oldPars = newPars
		newPars = c(parsGating, parsExperts)

print(summary(oldPars))
print(summary(newPars))

		if (sum((oldPars - newPars)^2) < tol)
			break

		## E step	
		if (J == 2) {
			predGating = predict(gating, newx = X[,colsGating, drop = FALSE], s = lambda[J+1], type = "response",
				offset = offsetGating)							# N x 1 matrix
			## transform to probs
			predGating = logistic(predGating)					# probabilities for group 2
			predGating = cbind(1 - predGating, predGating)		# N x 2 matrix
		} else {
			predGating = sapply(1:J, function(j) predict(gating[[j]], newx = X[,colsExperts, drop = FALSE], s = lambda[J+1],
				type = "response", offset = offsetExperts))		# N x J matrix
			## transform to probs
			predGating = softmax(predGating)					# N x J matrix (with probabilities for J groups)
		}
		if (K == 2) {
			predExperts = sapply(1:J, function(j) predict(experts[[j]], newx = X[,colsExperts, drop = FALSE], s = lambda[j],
				type = "response", offset = offsetExperts))		# N x J matrix (with J probabilities for the higher class level)
			## transform to probs
			predExperts = apply(predExperts, 2, logistic)
			## get probabilities of actual class
			predExperts[y == lev1[1],] = 1 - predExperts[y == lev1[1],]
		} else {	
			predExperts = list()
			## list of J N x K matrices
			for (j in 1:J) {
				predExperts[[j]] = sapply(1:K, function(k) predict(experts[[k]][[j]], newx = X[,colsExperts, drop = FALSE],
					s = lambda[j], type = "response", offset = offsetExperts))			# N x J matrix
			}
			## transform to probs
			predExperts = lapply(predExperts, softmax)
			## get probabilities of actual class
			predExperts = sapply(predExperts, function(pr) pr[cbind(1:N,as.numeric(y))])	# N x J matrix
		}
		weights = predGating * predExperts
		weights = weights/rowSums(weights)

print(summary(weights))
print(i)

		## transform weights
		weight.targets = weights
		weight.targets[weights < .Machine$double.eps] = .Machine$double.eps ## literature says 1e-03 works well
		## in the binary case keep only the second group
		weight.targets = log(weight.targets)
		if (J == 2) {
			weight.targets = weight.targets[,2]
		}
	}
	
	res = list(gating = gating, experts = experts, weights = weights, iters = i, J = J, K = K, colsGating = colsGating,
		colsExperts = colsExperts, lambda = lambda, offsetGating = offsetGating, offsetExperts = offsetExperts,
		lev1 = lev1, lev = lev)
	class(res) = "EMlnet"
	return(res)
}



#' @title Predict New Observations by a Trained Logistic Mixture of Experts Model
#'
#' @description
#' Predict new observations by a trained logistic mixture of experts model.
#'
#' @param object An object of class \code{EMlnet}.
#' @param newdata A \code{data.frame} or \code{matrix} with data for which to make predictions.
#' @param ... Further arguments. Currently unused.
#'
#' @return
#' A \code{list} with components:
#' \item{class}{A \code{factor} with predicted class levels.}
#' \item{posterior}{A \code{matrix} with predicted posterior probabilities.}
#' \item{gating}{The probabilities predicted by the fating model.}
#' \item{experts}{The class posterior probabilities predicted by individual experts.}
#'
#' @export

predict.EMlnet = function(object, newdata, ...) {
	if (!inherits(object, "EMlnet"))
		stop("'object' not of class 'EMlnet'")
	if (is.null(dim(newdata)))		# a vector is treated as one row
		dim(newdata) = c(1, length(newdata))
	X = as.matrix(newdata)
	if (object$J == 2) {
		predGating = predict(object$gating, newx = X[,object$colsGating, drop = FALSE], s = object$lambda[object$J+1],
			type = "response", offset = object$offsetGating)		# N x 1 matrix
		## transform to probs
		predGating = logistic(predGating)							# probabilities for group 2
		if (object$K == 2) {
			predExperts = sapply(1:object$J, function(j) predict(object$experts[[j]], newx = X[,object$colsExperts, drop = FALSE],
				s = object$lambda[j], type = "response", offset = object$offsetExperts))
				# N x J matrix (with J probabilities for the higher class level)
			## transform to probs
			predExperts = apply(predExperts, 2, logistic)	# N x 2 matrix
			## predGating: N x 1 matrix (group 2)
			## predExperts: N x 2 matrix (probs for class 2 from groups 1 and 2)
			post = rowSums(cbind(1 - predGating, predGating) * predExperts)			# N x 1 matrix (class 2)
			post = cbind(1 - post, post)											# N x 2 matrix (classes 1 and 2)
		} else if (object$K > 2) {
			predExperts = list()
			## list of J N x K matrices
			for (j in 1:object$J) {
				predExperts[[j]] = sapply(1:object$K, function(k) predict(object$experts[[k]][[j]],
					newx = X[,object$colsExperts, drop = FALSE], s = object$lambda[j], type = "response",
					offset = object$offsetExperts))	# N x J matrix
			}
			## transform to probs
			predExperts = lapply(predExperts, softmax)
			post = (1 - as.numeric(predGating)) * predExperts[[1]] + as.numeric(predGating) * predExperts[[2]]
		} else
			stop("something is wrong with 'object$K'")
	} else if (object$J > 2) {
		predGating = sapply(1:object$J, function(j) predict(object$gating[[j]], newx = X[,object$colsExperts, drop = FALSE],
			s = object$lambda[object$J+1], type = "response", offset = object$offsetExperts))	# N x J matrix
		## transform to probs
		predGating = softmax(predGating)					# N x J matrix (with probabilities for J groups)
		if (object$K == 2) {
			predExperts = sapply(1:object$J, function(j) predict(object$experts[[j]], newx = X[,object$colsExperts, drop = FALSE],
				s = object$lambda[j], type = "response", offset = object$offsetExperts))
			## transform to probs
			predExperts = apply(predExperts, 2, logistic)	# N x J matrix (with J probabilities for the higher class level)
			post = rowSums(predGating * predExperts)
			post = cbind(1 - post, post)
		} else if (object$K > 2) {
			predExperts = list()
			## list of J N x K matrices
			for (j in 1:object$J) {
				predExperts[[j]] = sapply(1:object$K, function(k) predict(object$experts[[k]][[j]],
					newx = X[,object$colsExperts, drop = FALSE], s = object$lambda[j], type = "response",
					offset = object$offsetExperts))			# N x J matrix
			}
			## transform to probs
			predExperts = lapply(predExperts, softmax)
			post = matrix(rowSums(sapply(1:object$J, function(z) predGating[,z] * predExperts[[z]])), ncol = object$K)
			# N x K matrix
		} else
			stop("something is wrong with 'object$K'")
	} else
		stop("Something is wrong with 'object$J'")
	rownames(post) = rownames(X)
	colnames(post) = object$lev1		
	cl = factor(object$lev1[max.col(post)], levels = object$lev)
	return(list(class = cl, posterior = post, gating = predGating, experts = predExperts))
}



#' @noRd
logistic = function(eta) {
	binomial()$linkinv(eta)
}



#' @noRd
softmax = function(eta) {
	maxeta = apply(eta, 1, max)
	eta = eta - maxeta
	probs = exp(eta)
	probs = probs/rowSums(probs)
}
