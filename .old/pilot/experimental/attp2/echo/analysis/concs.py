import uv

data = uv.BM3('20210303-concs.csv')
print(data.concs)
concs = data.concs / 0.005

def v1(c1, c2, v2):
    return (c2 * v2) / c1 

print(v1(concs, 10, 30_000)) # µM, µM, µl
