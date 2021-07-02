# round 2

## aim

measure effect of 0.1% bsa on:
- assay/signal stability over time @ room temp
- kd

## this time:
only selecting known binders to bm3 a82f/f87v. going to get those from lauras data and select 48 structurally diverse mutants.
I have about 700µl of bm3 a82f/f87v left over from round 1 (9th June 2021) at a conc of 1104 µM.

```python
In [2]: def v2(c1,v1,c2):
   ...:     return (c1 * v1) / c2
   ...:

In [3]: v2(1104, 700, 10)
Out[3]: 77280.0

In [5]: v2(1104, 700, 10) / (30 * 384)
Out[5]: 6.708333333333333
```

That's enough for around 6 full plates. This experiment would use two and the remaining 4 can be used for full library hit/miss screening.

```python
In [6]: def v1(c1,v2,c2):
   ...:     return (c2 * v2) / c1
   ...:

In [7]: v2 = 380*30 + 2

In [8]: v1(1104, v2, 10)
Out[8]: 103.27898550724638

In [9]: v2
Out[9]: 11402 # µl

In [10]: v2 / 1000
Out[10]: 11.402 # ml
```


- p1 : kpi-bm3
- p2 : bsa-bm3
- p3 : kpi
- p4 : bsa
