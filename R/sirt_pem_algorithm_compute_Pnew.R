## File Name: sirt_pem_algorithm_compute_Pnew.R
## File Version: 0.03

sirt_pem_algorithm_compute_Pnew <- function( tt, P0, P1, P2)
{
    Pnew <- (1-tt)^2 * P0 + 2*tt*(1-tt)*P1 + tt^2 * P2
    return(Pnew)
}
