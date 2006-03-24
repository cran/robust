stop.on.bdObject <- function(x)
{
	invisible()
}

key <- function(x, y, text, lines, corner = c(0, 1), cex = par("cex"),
	col = par("col"), lty = par("lty"), lwd = par("lwd"), ...)
{
	if(missing(x))
		x <- "topright"

	if(!is.null(lines$col))
		col = lines$col

	if(!is.null(lines$lty))
		lty = lines$lty

	if(!is.null(lines$lwd))
		lwd = lines$lwd

	legend(x = x, y = y, legend = unlist(text), xjust = corner[1], yjust = corner[2],
		col = col, lty = lty, lwd = lwd, bty = "n")

	invisible()
}
