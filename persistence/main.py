from filtrations import Filtration, ClassicFiltrations
from bar_code import BoundaryMatrix

filename = "filtrations/filtration_B.txt"

filtration = Filtration(filename)

bound_mat = BoundaryMatrix(filtration)

bound_mat.bar_codes(0.05)


# classic_filtrations = ClassicFiltrations([1, 2, 3, 10])
# classic_filtrations.filtrations_to_file("filtrations/classic/")