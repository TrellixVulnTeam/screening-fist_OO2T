# bsa stability testing
## aim
measure effect of 0.1% bsa buffer additive to:
- assay time stability 
- kd measurements for different compounds


## background
bsa can stibilize proteins in solution whilst being relatively inert, which is good for my assay because it has to be >20°C to avoid dmso freezing and causing scattering effects that mess with the measurements. I want to know how much it stabilizes the assay - how long does a plate last?

I also need to know if BSA binds to any of the compounds which can skew the readings

## experiments
### ```r1/``` 
First round - I picked 24 structurally diverse compounds from the fda library and screened them against BM3 a82F/F87V. I didn't check if they are already binders, which limits the usefulness of tis experiment. Also, a few of the compounds were low, so there will be some missing data. I can find out which ones by analyzing the exceptions report. I'll build a method for that into ```echo```.


### ```r2/```  
This time I will use Laura's data to select a set of compounds from known binders. I should try to do 48 compounds.
