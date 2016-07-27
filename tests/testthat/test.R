context("EMglmnet")

test_that("J = 2, K = 3", {
	X = as.matrix(iris[1:4])
	y = iris$Species
	lev3 = levels(y)
	
	res32 = EMglmnet(y, X, lambda = 1e-03)
	pred32 = predict(res32, X)

	expect_equal(levels(pred32$class), lev3)
	expect_equal(length(pred32$class), nrow(X))

	expect_equal(ncol(pred32$posterior), 3)
	expect_equal(colnames(pred32$posterior), lev3)
	expect_equal(nrow(pred32$posterior), nrow(X))
	expect_equal(rowSums(pred32$posterior), rep(1, nrow(X)))

	expect_equal(ncol(pred32$gating), 1)
	expect_equivalent(nrow(pred32$gating), nrow(X))
	
	for (j in 1:2) {
		expect_equal(ncol(pred32$experts[[j]]), 3)
		expect_equal(nrow(pred32$experts[[j]]), nrow(X))
		expect_equal(colnames(pred32$experts[[j]]), lev3)
	}
})	

test_that("J = 3, K = 3", {
	X = as.matrix(iris[1:4])
	y = iris$Species
	lev3 = levels(y)
	
	res33 = EMglmnet(y, X, lambda = 1e-03, J = 3)
	pred33 = predict(res33, X)
	
	expect_equal(levels(pred33$class), lev3)
	expect_equal(length(pred33$class), nrow(X))
	
	expect_equal(ncol(pred33$posterior), 3)
	expect_equal(colnames(pred33$posterior), lev3)
	expect_equal(nrow(pred33$posterior), nrow(X))
	expect_equal(rowSums(pred33$posterior), rep(1, nrow(X)))

	expect_equal(ncol(pred33$gating), 3)
	expect_equal(nrow(pred33$gating), nrow(X))
	expect_equivalent(rowSums(pred33$gating), rep(1, nrow(X)))

	for (j in 1:3) {
		expect_equal(ncol(pred33$experts[[j]]), 3)
		expect_equal(nrow(pred33$experts[[j]]), nrow(X))
		expect_equal(colnames(pred33$experts[[j]]), lev3)
	}
})	

test_that("J = 2, K = 2", {
	X = as.matrix(iris[1:100,1:4])
	y = iris$Species
	lev2 = levels(y)[1:2]
	y = factor(y[1:100], levels = lev2)

	res22 = EMglmnet(y, X, lambda = 1e-03)
	pred22 = predict(res22, X)
	
	expect_equal(levels(pred22$class), lev2)
	expect_equal(length(pred22$class), nrow(X))
	
	expect_equal(ncol(pred22$posterior), 2)
	expect_equal(colnames(pred22$posterior), lev2)
	expect_equal(nrow(pred22$posterior), nrow(X))
	expect_equivalent(rowSums(pred22$posterior), rep(1, nrow(X)))
	
	expect_equal(ncol(pred22$gating), 1)
	expect_equal(nrow(pred22$gating), nrow(X))
	
	expect_equal(ncol(pred22$experts), 2)	# only one class
	expect_equal(nrow(pred22$experts), nrow(X))
})	

test_that("J = 3, K = 2", {
	X = as.matrix(iris[1:100,1:4])
	y = iris$Species
	lev2 = levels(y)[1:2]
	y = factor(y[1:100], levels = lev2)

	res23 = EMglmnet(y, X, lambda = 1e-03, J = 3)
	pred23 = predict(res23, X)
	
	expect_equal(levels(pred23$class), lev2)
	expect_equal(length(pred23$class), nrow(X))
	
	expect_equal(ncol(pred23$posterior), 2)
	expect_equal(colnames(pred23$posterior), lev2)
	expect_equal(nrow(pred23$posterior), nrow(X))
	expect_equivalent(rowSums(pred23$posterior), rep(1, nrow(X)))
	
	expect_equal(ncol(pred23$gating), 3)
	expect_equal(nrow(pred23$gating), nrow(X))
	expect_equivalent(rowSums(pred23$gating), rep(1, nrow(X)))
	
	expect_equal(ncol(pred23$experts), 3)	# only one class
	expect_equal(nrow(pred23$experts), nrow(X))
})

## =================================================================================================================================
context("EMlnet")

test_that("J = 2, K = 3", {
	X = as.matrix(iris[1:4])
	y = iris$Species
	lev3 = levels(y)
	
	res32 = EMlnet(y, X, lambda = 1e-03)
	pred32 = predict(res32, X)

	expect_equal(levels(pred32$class), lev3)
	expect_equal(length(pred32$class), nrow(X))

	expect_equal(ncol(pred32$posterior), 3)
	expect_equal(colnames(pred32$posterior), lev3)
	expect_equal(nrow(pred32$posterior), nrow(X))
	expect_equal(rowSums(pred32$posterior), rep(1, nrow(X)))

	expect_equal(ncol(pred32$gating), 1)
	expect_equivalent(nrow(pred32$gating), nrow(X))
	
	for (j in 1:2) {
		expect_equal(ncol(pred32$experts[[j]]), 3)
		expect_equal(nrow(pred32$experts[[j]]), nrow(X))
		# expect_equal(colnames(pred32$experts[[j]]), lev3)
	}
})	

test_that("J = 3, K = 3", {
	X = as.matrix(iris[1:4])
	y = iris$Species
	lev3 = levels(y)
	
	res33 = EMlnet(y, X, lambda = 1e-03, J = 3)
	pred33 = predict(res33, X)
	
	expect_equal(levels(pred33$class), lev3)
	expect_equal(length(pred33$class), nrow(X))
	
	expect_equal(ncol(pred33$posterior), 3)
	expect_equal(colnames(pred33$posterior), lev3)
	expect_equal(nrow(pred33$posterior), nrow(X))
	expect_equal(rowSums(pred33$posterior), rep(1, nrow(X)))

	expect_equal(ncol(pred33$gating), 3)
	expect_equal(nrow(pred33$gating), nrow(X))
	expect_equivalent(rowSums(pred33$gating), rep(1, nrow(X)))

	for (j in 1:3) {
		expect_equal(ncol(pred33$experts[[j]]), 3)
		expect_equal(nrow(pred33$experts[[j]]), nrow(X))
		# expect_equal(colnames(pred33$experts[[j]]), lev3)
	}
})	

test_that("J = 2, K = 2", {
	X = as.matrix(iris[1:100,1:4])
	y = iris$Species
	lev2 = levels(y)[1:2]
	y = factor(y[1:100], levels = lev2)

	res22 = EMlnet(y, X, lambda = 1e-03)
	pred22 = predict(res22, X)
	
	expect_equal(levels(pred22$class), lev2)
	expect_equal(length(pred22$class), nrow(X))
	
	expect_equal(ncol(pred22$posterior), 2)
	expect_equal(colnames(pred22$posterior), lev2)
	expect_equal(nrow(pred22$posterior), nrow(X))
	expect_equivalent(rowSums(pred22$posterior), rep(1, nrow(X)))
	
	expect_equal(ncol(pred22$gating), 1)
	expect_equal(nrow(pred22$gating), nrow(X))
	
	expect_equal(ncol(pred22$experts), 2)	# only one class
	expect_equal(nrow(pred22$experts), nrow(X))
})	

test_that("J = 3, K = 2", {
	X = as.matrix(iris[1:100,1:4])
	y = iris$Species
	lev2 = levels(y)[1:2]
	y = factor(y[1:100], levels = lev2)

	res23 = EMlnet(y, X, lambda = 1e-03, J = 3)
	pred23 = predict(res23, X)
	
	expect_equal(levels(pred23$class), lev2)
	expect_equal(length(pred23$class), nrow(X))
	
	expect_equal(ncol(pred23$posterior), 2)
	expect_equal(colnames(pred23$posterior), lev2)
	expect_equal(nrow(pred23$posterior), nrow(X))
	expect_equivalent(rowSums(pred23$posterior), rep(1, nrow(X)))
	
	expect_equal(ncol(pred23$gating), 3)
	expect_equal(nrow(pred23$gating), nrow(X))
	expect_equivalent(rowSums(pred23$gating), rep(1, nrow(X)))
	
	expect_equal(ncol(pred23$experts), 3)	# only one class
	expect_equal(nrow(pred23$experts), nrow(X))
})
