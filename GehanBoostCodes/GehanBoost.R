ngrad.Gehan <- function(y,f,w=1) {
	U <- as.vector(y[,1])
	Del <- as.vector(y[,2])
	f <- as.vector(f)
	res <- U - f
	n <- length(U)
	tmpg.fun <- function(I,res,del) {
		ei <- res[I]
		diff <- ei - res
		arg2 <- sum(del*(diff>=0))
		
		if(del[I] > 0.5) arg1 <- sum(diff<=0)
		else arg1 <- 0
		
		arg1 - arg2
		}
	g <- unlist(lapply(1:n,tmpg.fun,res,Del))
	-1*g/n 
	}
rho.Gehan <- function(y,f,w=1) {
	U <- as.vector(y[,1])
	Del <- as.vector(y[,2])
	f <- as.vector(f)
	res <- U - f
	n <- length(U)
	tmp.fun <- function(I,res,del) {
		if(del[I] < 0.5) out <- 0
		else {
			ei <- res[I]
			diff <- ei - res
			out <- -1 * sum(diff[diff <= 0])
			}
		out
		}
	out <- unlist(lapply(1:n,tmp.fun,res,Del))/n
	out
	}
Gehan <- Family(ngradient=ngrad.Gehan,loss=rho.Gehan,
	check_y = function(y) {
		if (!inherits(y,"Surv")) stop("response is not an object of class ", sQuote("Surv"), " but ", sQuote("family = Gehan()"))
		TRUE
		},
	name = "Gehan loss")


