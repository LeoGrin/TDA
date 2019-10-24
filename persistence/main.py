from filtrations import Filtration, ClassicFiltrations
from bar_code import BoundaryMatrix
import matplotlib.pyplot as plt
import numpy as np


#compute filtrations for classic topological spaces
classic_filtrations = ClassicFiltrations([1, 2, 3, 10])
classic_filtrations.filtrations_to_file("filtrations/classic/")

l_p = []
l_t = []

names = ["filtration_B", "filtration_C", "filtration_A", "filtration_D"]

for name in classic_filtrations.all_names:
    print("computing " + name)
    filtration = Filtration("filtrations/classic/" + name + ".txt")
    bound_mat = BoundaryMatrix(filtration)
    res = bound_mat.bar_codes(0)
    l_t.append(res[1])
    l_p.append(res[2])

for name in names:
    print("computing " + name)
    filtration = Filtration("filtrations/" + name + ".txt")
    bound_mat = BoundaryMatrix(filtration)
    if name is "filtration_B":
        res = bound_mat.bar_codes(0.05)
    else:
        res = bound_mat.bar_codes(0)
    l_t.append(res[1])
    l_p.append(res[2])





x = np.log(l_p)
y = np.log(l_t)

deg = 1
coef = np.polyfit(x, y, deg=deg)[::-1]
print(coef)
y_poly = [sum([coef[i] * elem ** i for i in range(deg+1)]) for elem in x]

plt.scatter(x, y, marker="o", c="b")
plt.plot(x, y_poly, c="r")
plt.xlabel("log(n_simplex)")
plt.ylabel("log(time)")
plt.legend("pente = {}".format(coef[1]))
plt.show()