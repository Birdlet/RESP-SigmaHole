# convert psi4 hessian to xtb hessian format
import sys

inf = sys.argv[1]
outf = sys.argv[2]

xss = []
with open(inf) as f:
    lines = list(f.readlines())
    cs = lines[0].strip().split()
    nlines = int(cs[0]) * int(cs[1])
    for l in lines[1:1+nlines]:
        xs = [float(x) for x in l.strip().split()]
        xss.extend(xs)

natoms = int(cs[0])
nhess = int(cs[1])

with open(outf, 'w') as f:
    f.write("$hessian")
    for i in range(nhess):
        for j, x in enumerate(xss[i*nhess:(i+1)*nhess]):
            if j % 5 == 0:
                f.write("\n%20.10f" % x)
            else:
                f.write("%15.10f" % x)

"""

set hessian_write on

H, wfn = hessian('b3lyp-d3', return_wfn=True)

wfn.hessian().print_out()

np.array(H)
# cubeprop(wfn)
"""
