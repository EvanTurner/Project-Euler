#! /usr/bin/env python
# Project Euler problem 101: Optimum polynomial

def u(n, m): return sum([(-1)**i * n**i for i in xrange(m+1)])

# takes n + 1 data points -- (x, f(x)) -- and generates interpolated poly of order n
def get_newton_diffs(x, fx):
    n = len(x) - 1

    diffs = []
    
    diff_list = fx              # diff level 0
    diffs.append(diff_list)     # add this set of divided differences
    diff_list_prev = diff_list
    
    diff_level = 1
    
    while (diff_level <= n):
        diff_list = []
        i, j = 1, 0
        while i <= n - (diff_level-1):
            diff_list.append( 
                    (diff_list_prev[i] - diff_list_prev[j]) / (x[i + (diff_level-1)] - x[j])
                )
            i += 1
            j += 1            
        
        diffs.append(diff_list)     # add this set of divided differences
        diff_list_prev = diff_list
        diff_level += 1

    return diffs

def eval_newton_poly(diffs, x_vals, eval=0):
    n = len(x_vals) - 1  # n + 1 data points for polynomial order n
    fx_eval = 0
    
    for i in xrange(0, n+1):
        coeff = 1
        for j in xrange(i):
            coeff *= (eval - x_vals[j])
        fx_eval += coeff * diffs[i][0]
    
    return fx_eval

def generating_fncn(x, fncn=0):
    f_x = 0

    if fncn==0:
        f_x = u(x, 10)  # pr101 default generating function
    elif fncn==1:
        f_x = x ** 3    # cubic function, for testing

    return f_x

def OP(k,n, diffs): return eval_newton_poly(diffs, k, n)

def solve_pr101(k_max, fncn=0):
    epsilon = 1e-03
    total_FIT = 0
    
    for i in xrange(1, k_max+1):
        k = xrange(1, i+1)
        f_k = [generating_fncn(x, fncn) for x in k]
        diffs = get_newton_diffs(k, f_k)
        
        for x in xrange(len(k), len(k) + 20):
            eval_newton_x = eval_newton_poly(diffs, k, x)
            gf_x = generating_fncn(x, fncn)
            if abs(eval_newton_x  - gf_x) > epsilon:
                total_FIT += eval_newton_x
                break


    return total_FIT
