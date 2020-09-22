.ode_compute_u <- function(alpha, beta, u_0, tau){
    u_0*exp(-beta*tau) + alpha/beta*(1-exp(-beta*tau))
}

.ode_compute_s <- function(alpha, beta, gamma, u_0, s_0, tau){
    s_0*exp(-gamma*tau) + alpha/gamma*(1-exp(-gamma*tau)) + (alpha - beta*u_0)/(gamma - beta)*(exp(-gamma*tau) - exp(-beta*tau))
}

.ode_simulate_discrete <- function(param_df, u_0, s_0){
    stopifnot(ncol(param_df) == 4, all(colnames(param_df) %in% c("t", "alpha", "beta", "gamma")))
    len <- nrow(param_df)
    
    res_df <- data.frame(t = param_df$t, u = rep(NA, len), s = rep(NA, len))
    res_df[1,"u"] <- u_0
    res_df[1,"s"] <- s_0
    
    for(i in 2:len){
        alpha <- param_df$alpha[i-1]
        beta <- param_df$beta[i-1]
        gamma <- param_df$gamma[i-1]
        t_diff <- res_df$t[i]-res_df$t[i-1]
        u_prev <- res_df$u[i-1]
        s_prev <- res_df$s[i-1]
        
        res_df$u[i] <- u_prev + (alpha - beta*u_prev)*t_diff
        res_df$s[i] <- s_prev + (beta*u_prev - gamma*s_prev)*t_diff
    }
    
    res_df
}

.ode_simulate_closed <- function(param_df, u_0, s_0){
    stopifnot(ncol(param_df) == 4, all(colnames(param_df) %in% c("t", "alpha", "beta", "gamma")))
    len <- nrow(param_df)
    
    res_df <- data.frame(t = param_df$t, u = rep(NA, len), s = rep(NA, len))
    res_df[1,"u"] <- u_0
    res_df[1,"s"] <- s_0
    
    changepoints <- c(1,which(abs(diff(param_df$alpha)) >= 1e-6), len)
    
    for(phase in 1:(length(changepoints)-1)){
        t_0 <- res_df$t[changepoints[phase]]
        u_0 <- res_df$u[changepoints[phase]]
        s_0 <- res_df$s[changepoints[phase]]
        
        for(i in (changepoints[phase]+1):changepoints[phase+1]){
            alpha <- param_df$alpha[i-1]
            beta <- param_df$beta[i-1]
            gamma <- param_df$gamma[i-1]
            tau <- param_df$t[i] - t_0
            
            res_df$u[i] <- .ode_compute_u(alpha, beta, u_0, tau)
            res_df$s[i] <- .ode_compute_s(alpha, beta, gamma, u_0, s_0, tau)
        }
    }
    
    res_df
}

.ode_steadystate_time <- function(alpha, beta, gamma, u_0, s_0, epsilon = 0.03, 
                                  t_min = 0.01, t_max = 1000){
    # determine (u,s) at various times
    t_seq <- exp(seq(log(t_min), log(t_max)), length.out = 20)
    position_mat <- sapply(t_seq, function(t_val){
        
        names(vec) <- c("u", "s")
        vec
    })
    
    eval_func <- function(main_arg, alpha, beta, gamma, epsilon){
        u <- .ode_compute_u(alpha, beta, u_0, t_val)
        s <- .ode_compute_u(alpha, beta, gamma, u_0, s_0, t_val)
        abs(beta/gamma*u - s) <= epsilon
    }
    
    # determine which ones are close to steady state
    bool_vec <- sapply(t_seq, function(t_val){
        eval_func(main_arg = t_val, alpha, beta, gamma, epsilon)
    })
    idx <- length(bool_vec)
    stopifnot(bool_vec[idx])
    while(bool_vec[idx-1]){
        idx <- idx-1
    }
    t_high <- t_seq[idx]
    t_low <- t_seq[idx-1]
    
    # binary search
}

.binary_search <- function(eval_func, hi, lo, max_iter = 10, ...){
    iter <- 1
    hi_bool <- eval_func(main_arg = hi, ...)
    lo_bool <- eval_func(main_arg = lo, ...)
    stopifnot(hi_bool != lo_bool)
    
    while(iter <= max_iter){
        mid <- (hi+lo)/2
        mid_bool <- eval_func(main_arg = mid, ...)
        
        if(mid_bool == hi_bool){
            hi <- mid
        } else {
            lo <- mid
        }
    }
    
    mid
}