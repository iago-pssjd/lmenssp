variogram <-
function(resid, timeVar, id, binwidth, numElems = 0){

v  <- unlist(lapply (tapply (resid, id, function (x) outer (x, x, function (x, y) (x - y)^2/2)), function (x) x[lower.tri (x)]))
u  <- unlist(lapply (tapply (timeVar, id, function (x) outer (x, x, function (x, y) abs(x - y))), function (x) x[lower.tri (x)]))

bins <- seq(floor(min(u)), ceiling(max(u)), binwidth)

bin.means <- NULL
bin.sizes <- NULL

for (i in 1 : (length(bins) - 1)){

bin.means <- c(bin.means, mean(cbind(u, v)[cbind(u, v)[, 1] >= bins[i] & cbind(u, v)[, 1] < bins[i+1], 2]))
bin.sizes <- c(bin.sizes, length(cbind(u, v)[cbind(u, v)[, 1] >= bins[i] & cbind(u, v)[, 1] < bins[i+1], 2]))

}

bin.mids <- NULL
for (i in 1 : (length(bins) - 1)) bin.mids <- c(bin.mids, (bins[i] + bins[i+1]) / 2)

plot(bin.mids[which(bin.sizes > numElems)], bin.means[which(bin.sizes > numElems)], 
     type = "l", xlab = "Lag", ylab = "Semivariogram")

output           <- list()
output$bin.mids  <- bin.mids
output$bin.means <- bin.means
output$bin.sizes <- bin.sizes

output

}
