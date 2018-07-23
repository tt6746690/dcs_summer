


```sh
matlab -nodesktop
```



```
S=5
>>>>> fminsearch (un-constrained):

mu =

   -0.1377
   -0.9931
   -2.4841
    1.3091
    1.9612

>>>>> fmincon (constrained):

mu =

    0.6243
   -0.4541
   -2.2364
    1.7449
    1.9662

hessian = [
    1.8779   -0.7298   -5.3524    3.8563    4.2885    0.0160
   -0.7298    0.5267    2.4054   -1.6905   -1.9628   -0.0012
   -5.3524    2.4054   17.8072  -12.4082  -13.6724   -0.0518
    3.8563   -1.6905  -12.4082    9.0886    9.7160    0.0356
    4.2885   -1.9628  -13.6724    9.7160   11.2498    0.0308
    0.0160   -0.0012   -0.0518    0.0356    0.0308    0.0004
]
det(hessian) = 3.2425e-06 > 0 OK!

>>>>> fmincon with hessian (constrained):

mu =

    1.7370
   -0.8750
   -2.2091
   -0.7607
   -1.8342


>>>>> gradient ascent (constrained):

mu =

    0.0250
   -2.0751
   -2.6266
    0.6562
    1.3252
```




S=1
```
>>>>> fminsearch (un-constrained):

mu =

   -2.0650

>>>>> fmincon (constrained):

mu =

    1.9987

>>>>> fmincon with hessian (constrained):

mu =

   -1.2710

>>>>> gradient ascent (constrained):

mu =

    2.0471
```



S=2
```
>>>>> fmincon (constrained):

mu =

   -1.7024
    1.8452
```




