

si_fmt <- function(x) { 
	system2(here( "bin/si-format.sh"), x, stdout = TRUE) 
	}