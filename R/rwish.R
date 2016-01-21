rwish <-
function (S, df) {
	.Call( "rwish", S, df:(df - ncol(S) + 1), PACKAGE = "BayesComm" )
}
