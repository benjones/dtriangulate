
# Predicate Derivations

Analysis of the error bounds for predicates making use of FMA, using the techniques from Shewchuk


## Cancellation and FMA

```
R = (ax - cx)(by - cy) - (ay - cy)(bx - cx)
expand

R = axby - axcy - bycx + cxcy - aybx + aycx + bxcy - cycx //cxcy terms cancel

R = ax(by - cy) - ay(bx - cx) + bxcy - bycx

R = fma(ax, by - cy, bxcy) - fma(ay, bx - cx, bycx) 

```

For analysis, lowercase == true value, uppercase == FP approximations, (?) == FP version of ?:

```
d = by - cy; D = by (-) cy ; d = D +- eps|D|
e = cx - bx; ...
f = bycx; ...
g = bxcy; ...

h = ax*(by - cy) + bxcy = ax*d + g;  H = (ax * D + G)( 1 +- eps) //one round for FMA
ax*D + G = H +- eps|H|
//substitute D, G:

h = ax(D +/ eps|D|) + G +/ eps|G|

h = ax*D + G +-eps(|D| + |G|)
//sub H for the first term

h = H +- eps(|H| + |D| + |G|)


similarly,
k = fma(ay, bx - cy, bycx) = K +- eps(|K| + |E| + |F|)
```

Now combine h and k:

```
det = h - k;  DET = (H - K)(1 +- eps);  H - K = DET(1 +- eps)
det = H - K +-eps(|D| + |E| + |F| + |G| + |H| + |K|)

det = DET +- eps|DET| +- eps(|D| + |E| + |F| + |G| + |H| + |K|)

```

The sign is correct if `(1 - eps)|Det| > eps(|D| + |E| + |F| + |G| + |H| + |K|)`

Dividing by `(1 - eps)` gives

```
|DET| >= (eps + eps^2 + 2eps^3)(|D| + |E| + |F| + |G| + |H| + |K|) 

```

The epsilon term (I think) can be expanded as far as you want, and for some power p, you can just bound it by 2eps^p, but I just stopped at eps^3 since even it will underflow when added to eps

Let `s` be the sum of absolute values:

`s = S +- 6eps|S|` where 6 is the number of terms

So we know the sign if 

```
|DET| >= (eps + eps^2 + 2eps^3)(S + 6 eps S)
|DET| >= (7eps + 7eps^2 + 9eps^3)S
|DET| >= (7eps + 8eps^2)S
```
This seems to be worse than using the factored representation directly: `|A| >= (3 eps + 16eps^2)(|acx*bcy| + |acy*bcx|)`
Especially since '|S| >= ' the terms summed the Shewchuk version.

Testing on random inputs, this version gets the wrong answer more often than the naive code (first definition of R above).  
That said, both were wrong on the order of 1 in 1 million random sets of points.


## Just FMA the last step

```
R = FMA(ax - cx, by - cy, (cy - ay)(bx - cx));

d = ax - cy
e = by - cy
f = cy - ay
g = bx - cx

D = ax(-) cy   -> d = D(1 +- eps)
Same for E, F, G

h = fg; H = G (*) G = FG(1 +- eps)
H = FG +- eps|H|

h = (F +- eps|F|)(G +- eps|G|)
h = FG +-(2eps + eps^2)|FG|
h = H +- eps|H| +- (2eps + eps^2)|FG|
h = H +- eps|H| +- (2eps + eps^2)(|H| +_ eps|H|)
h = H +- |H|(3eps +- 3eps^2 + eps^3)

det = de + h
DET = (DE + H)(1 +- eps)
det = (D +-eps|D|)(E _+-eps|E|) + H +- |H|(3eps +- 3eps^2 + eps^3)
det = DE + H + 2eps|DE| + eps^2|DE| + |H|(3eps +- 3eps^2 + eps^3)
det = DET +- eps|DET| + (2eps + eps^2)|DE| + |H|(3eps +- 3eps^2 + eps^3)

```

We don't have to compute DE by itself since it's part of the FMA, so we'd like to avoid including it in the error analysis:

```
det = DET +- eps|DET| + (2eps + eps^2)(|DE| + |H|) + (eps + 2eps^2 + eps^3)|H|

|DE| + |H| <= |DET| + 2|H| 

det = DET +- eps|DET| +- (2eps + eps^2)(|DET| + 2|H|) + (eps + 2eps^2 + eps^3)|H|

det = DET +- (3eps + eps^2)|DET| +- (5eps + 4eps^2 + eps^3)|H|

```

det has the right sign if 

```
(1 - 3eps - eps^2)|DET| > (5eps + 4eps^2 + eps^3)|H|
```

dividing by `(1 - 3eps - eps^2) on both sides gives

```
|DET| > (5eps + 19eps^2 + 63eps^3 + (208eps^4 + 63eps^5)/(1 - 3eps + eps^2) )|H|

|DET| > (5eps + 19eps^2 + 64eps^3)|H|
```

The leading term (5eps vs 3eps for the shewchuck version) is higher, but we're multiplying it by |H|, not |fl((ax - cy)(by - cy))| (+) |fl((ay - cy)(bx - cx))|.  
|H| is the second half only of that expression, so will be less than what Shewchuck is multiplying.  

In practice, I haven't seen a case where this version is incorrect, but the naive version is correct (I think it can be rigorously proven, but I haven't)

Since D, E, F, and G will be computed, we could compute bothe FMA(D, E, FG) and (F, G, DE).  
The error bound will be a function of either FG or DE, so pick the one with smaller absolute value.
