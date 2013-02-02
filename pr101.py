#! /usr/bin/env python
# Project Euler problem 101: Optimum polynomial

def u(n, m): return sum([(-1)**i * n**i for i in xrange(m+1)])

# takes n + 1 data points -- (x, f(x)) -- and generates divided differences
def get_newton_diffs(data_points):
    n = len(data_points) - 1    # degree of polynomial

    diffs = [[dp[1] for dp in data_points]]         # diff level 0
    
    for diff_level in xrange(1, n+1):
    
        diffs.append([])
        i, j = 1, 0

        while i <= n - (diff_level-1):
            diffs[diff_level].append(
                    (diffs[diff_level-1][i] - diffs[diff_level-1][j]) / 
                    (data_points[(diff_level-1) + i][0] - data_points[j][0])
                )
            i += 1
            j += 1     
        
    return diffs

def eval_newton_poly(diffs, data_points, eval=0):
    n = len(data_points) - 1  # n + 1 data points for polynomial order n
    fx_eval = 0
    
    for i in xrange(0, n+1):
        coeff = 1
        for j in xrange(i):
            coeff *= (eval - data_points[j][0])
        fx_eval += coeff * diffs[i][0]
    
    return fx_eval

def generating_fncn(x, fncn=0):
    f_x = 0

    if fncn==0:
        f_x = u(x, 10)  # pr101 default generating function
    elif fncn==1:
        f_x = x ** 3    # cubic function, for testing

    return f_x

def solve_pr101(k_max, fncn=0):
    epsilon = 1e-03
    max_n = 20
    total_FIT = 0
    
    for i in xrange(1, k_max+1):
    
        data_points = [(x, generating_fncn(x, fncn)) for x in xrange(1, i+1)]
        
        diffs = get_newton_diffs(data_points)
        n = len(data_points)
        
        for x in xrange(n, max_n):
            OP_k_n = eval_newton_poly(diffs, data_points, x)
            f_x = generating_fncn(x, fncn)
            
            if abs(OP_k_n  - f_x) > epsilon:
                total_FIT += OP_k_n   # This is the First Incorrect Term
                break

    return total_FIT
