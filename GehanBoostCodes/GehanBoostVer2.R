#
# Gehan family for 'mboost' package
#
# Brent Johnson
# 6 December 2011
#
# Reference : Johnson BA and Long Q (2011) Survival ensembles by the sum of pairwise differences with application to 
#		lung cancer microarray studies, Annals of Applied Statistics, 5:1081-1101.
#
# History : 
# 2009 : original Gehan boost function sets w=1 throughout
# 6 Dec 2011 : v2. rewrite for w = vector of integers, i.e.
#		"individual data is replicated as many time as indicated by weights"
#
# 7 Dec 2011 : 
#	getting R warning from cvrisk : In y * weights :
#  longer object length is not a multiple of shorter object length
#
# 9 Dec 2011 : 
#	note: default 'nu' in boost_control is often too small
#
Gehan <- Family(ngradient=function(y,f,w) {
	U <- as.vector(y[,1])
	Del <- as.vector(y[,2])
	n <- length(U)
	if(length(f)==1) f <- rep(f,n)
	if(length(w)==1) w <- rep(w,n)
	
	idx <- rep(1:n,w)
	Del <- Del[idx]
	n <- length(Del)
	res <- U[idx] - f[idx]
	tmpg.fun <- function(I,res,del) {
		ediff <- res[I] - res
		arg2 <- sum(del*(ediff>=0))
		arg1 <- ifelse(del[I] > 0.5,sum(ediff<=0),0)
		arg1 - arg2
		}
	g <- unlist(lapply(1:n,tmpg.fun,res,Del))
	-1*g/n},
	loss = gehan.loss <- function(y,f,w) {
		U <- as.vector(y[,1])
		Del <- as.vector(y[,2])
		n <- length(U)
		if(length(f)==1) rep(f,n)
		if(length(w)==1) rep(w,n)
		idx <- rep(1:n,w)		
		Del <- Del[idx]
		res <- U[idx] - f[idx]
		n <- length(Del)
		tmpl.fun <- function(I,res,del) {
			if(del[I] < 0.5) out <- 0
			else {
				ediff <- res[I] - res
				out <- -1 * sum(ediff[ediff <= 0])
				}
			out
			}
		out <- unlist(lapply(1:n,tmpl.fun,res,Del))/n
		out
		},
	risk = risk <- function(y, f, w = 1) sum(gehan.loss(y, f, w),na.rm=TRUE),
	offset = function(y,w) optimize(risk, interval = c(0,max(y[,1],na.rm=TRUE)),y=y,w=w)$minimum,
	check_y = function(y) {
		if (!inherits(y,"Surv")) stop("response is not an object of class ", sQuote("Surv"), " but ", sQuote("family = Gehan()"))
		y
		}, 
	name = "Gehan loss")
