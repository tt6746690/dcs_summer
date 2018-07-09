

##### 5.1

+ research paper
    + idea
        + probability of losses given a credit portfolio L
        + L: set of bonds, with interest
        + L is probabilistic model 
    + L(Z, \Sigma): 2 normal(0, 1) distribution 
        + first 2 chapters
        + biggest loss: everything defaults
        + also loss: reduce credit rating (AA -> A)
        + gain: increase credit rating
    + try to estimate with monte carlo 
    ```
    P(L(Z, \Sigma) > e) = E(1_{L(Z, \Sigma) > \epsilon})
    \approx \frac{1}{MN} \sum_{n=1}^N \sum_{n=1}^M 1_{L(Z, \Sigma) > \epsilon}
    ```
    + problem, MC runs too slowly
+ importance sampling 
    + how fast MC runs depends on variance
        + i..e Var_f(g(x))
    ```
    E(f(x)) = \int g(x) f(x) dx
    = \int \frac{g(x) f(x)}{\hat{f}(x)} \hat{f}(x) dx
    = E_{\hat{f}}{ \frac{g(x) f(x)}{\hat{f}(x)} }
    ```
    + potential source of error 
        + 2 MC on lhs and rhs, they should be approximately equal
        + but there seems to be some consistent bias, not as close as wed expect
    + importance sampling reduce variance 
        + var_{\hat{f}}{\frac{g(x) f(x)}{\hat{f}(x)}}
+ variance reduction method
    + why variance reduction for glasserman li doesnt work
+ 2-state

##### 5.21

+ stock return
    + log normal model
    + not that accurate, underestimate large losses
    + want thicker tails
+ motivation for H
    + P(y \leq H_{c(n)}^n)
+ why not model the covariance matrix and in this way sampling uses only 1 layer
    + 2 level useful
+ numerical analysis computing resources 
    + loginto appS0, appS1, ... then login to
    + compS0, compS1, ...
        + `ssh kjr@compS0.cs.toronto.edu`
+ TODOS:
    + check 2 simple MC algorithm 2.1 2.2 yield same result 
        + the result are consistently different, but try 
            + run on larger sample sizes 
            + try run 2 algorithm on the same saved parameters, maybe they are different and somehow, however unlikely, effect the result
        + also some MC estimate in 1 run yield 0 ??!!

##### 5.28

+ TODO:
    + verify 2.1 2.2 again with CI, on large sample sizes
    + re-read glasserman-li
+ since model is simulated, can compute the real statistics not relying on sampling
    + try different sample size -> observe varied variance ...


##### 6.10 

+ Readings 
    + glasserman and li 
        + likelihood function 
        + MC simulation 
    + chapter 3, (maybe 4)
+ Look at Glasserman Li impl
    + see if the likehood is adjusted ... properly in (Glasserman-LI - 21)


##### 6.25

+ look carefully at likelihood for the outer level change of normal distribution
    + see if expected value of likelihood sum up to 1
+ points 
    + changing `weights = EAD .* LGC` to `weights = LGC` only is not the solution 
        + varying in tail `l=0` to `l=0.99` doesn't change the shifted mean `mu`
        + but the original algorithm seems to have extreme `mu` (i.e. -2) even when `l=0`, (should expect no shift)
            + so a problem with finding the `mu`
            + TODO: test for this, using original weights, vary tail, expect change in shifted mu 
    + values of `mu`, i.e. shifted mean for `Z ~ N(mu, I_S)` effect MC estimates
        + maybe because number of samples for Z is not sufficient to offset the shift
        + so need to increase the number of samples for Z, does it fix the shifting problem?
    + glasserman&li works when not shifting `mu ~ 0`
        + corresponds to adam sturge's experiment

+ 2 sources of problem 
    + algorithm for computing shifted mean `mu`
        + effect variance reduction
    + likelihood value not accurate
        + regardless of values of `mu`, likelihood should compensate for the shift,
        + hypothesis: did not do sufficient sampling of `Z`, so that the likelihood does not compensate for the shift
            + yes, more sampling 

##### 7.2

1. even with the wrong `mu`, MC of likelihood should still converge to 1, for large enough sample size for `E` and `Z`
    + test with running the algorithm with large enough iterations, ... after setting the `mu` to 1, 2
2. Investigate more into how `mu` is picked,
    + potential problems
        + `mu` is not a local max, so the method didnt work 
            + verify this by conmputing the gradient and see if the 1st order derivative is 0 and 2nd order gradient <0
        + `mu` is a local max, but not the one we want
            + dunno how
    + run the algorithm with 
        + constarinted optimization 
        + unconstrainted steepest ascent
    + do some plotting in 1d/2d case

##### 7.4 meeting with kjr on presentation

+ presentation structure
    + MC + IS
        + introduction, definitions 
    + fiannce background
        + definitions
        + copula factor model
    + 2 lvl IS applied to copula factor model
        + rare event simulation ...
        + emphasis on why/how 2 lvl IS works 
        + sketch of algorithm 
    + research problem and approaches
        + downward bias of estimates with Glasserman&Li algorithm
        + need large enough sample for computing `mu` ?
        + `mu` did not find the right one ?
