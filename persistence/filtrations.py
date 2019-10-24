import itertools as itertools


class Simplex:
    def __init__(self, val, dim, vert):
        self.val = val
        self.dim = dim
        self.vert = vert

    def __eq__(self, simplex2):
        return self.vetices.__eq__(simplex2.verticles)

    def toString(self):
        return "val = {}; dim = {}; vert = {}".format(self.val, self.dim, "/".join(self.vert))


class Filtration:
    def __init__(self, filename):
        simplex_list = []
        with open(filename, "r") as test_filtration:
            for line in test_filtration.readlines():
                inputs = line.strip().split(" ")
                f, dim = float(inputs[0]), int(inputs[1])
                vertices = sorted(
                    list(inputs[2:]))  # vertices are sorted so we can compare between list of vertices
                simplex_list.append(Simplex(f, dim, vertices))

        simplex_list = sorted(simplex_list, key=lambda x: (
        x.val, x.dim))  # sort according to val and then to dim (to have a simplicial complex)

        dic_simplex = {}  # maps a simplex (as a sorted tuple of vertices) to its indentifier
        for i, s in enumerate(simplex_list):
            dic_simplex[tuple(s.vert)] = i

        self.dic_simplex = dic_simplex
        self.simplex_list = simplex_list


class ClassicFiltrations:
    """class to provide filtration of classic topological spaces"""
    def __init__(self, dimensions):
        """dimensions argument is for computing the filtrations of the ball and the sphere of given dimension"""
        self.all_filtrations = []
        self.all_names = []
        self.compute_filtrations(dimensions)

    def top_simplices_from_figure(self, list_vertices):
        # Useful function to list the top-dimensional subsimplices from the Moebius band, the torus and the Klein bottle.
        """
        Input: 2D-list from the vertices, row by row
        Output: set from the top-dimensional simplices (here the 2D simplices), i.e. set of list of vertices
        """
        top_simplices = []
        for i in range(len(list_vertices) - 1):
            for j in range(len(list_vertices[i]) - 1):
                top_simplices.append((list_vertices[i][j], list_vertices[i][j + 1], list_vertices[i + 1][j]))
                top_simplices.append((list_vertices[i][j + 1], list_vertices[i + 1][j], list_vertices[i + 1][j + 1]))
        return set(top_simplices)

    def subsimplices_from_simplex(self, simplex_vert):
        """
        Input: tuple from the vertices of the simplex
        Output: set from the subsimplices, i.e. set of tuple of vertices. Includes the input simplex.
        """
        return set(
            itertools.chain(*[itertools.combinations(simplex_vert, dim) for dim in range(1, len(simplex_vert) + 1)]))

    def all_simplices_from_top_simplices(self, simplices_set):
        """
        Input: set from the simplices, i.e. set of tuples of vertices
        Output: list from the simplices sorted in crescent order for the dimension and for the name
        """
        all_simplices = simplices_set.copy()
        for simplex_vert in simplices_set:
            new_simplices = self.subsimplices_from_simplex(simplex_vert)
            for new_simplex in new_simplices:
                all_simplices.add(new_simplex)
        all_simplices = sorted(list(all_simplices), key=lambda x: (len(x), x[0]))
        return all_simplices

    def from_simplices_to_filtration(self, sorted_simplices_list):
        """
        Input: list from the simplices sorted in crescent order for the dimension and for the name
        Output: filtration = list from the simplices sorted in crescent order for the dimension
        (according to the syntax
            Time of appearance in the filtration | Dimension | List from the vertices
        )
        In our filtration :
            Time of appearance(simplex K) = dim(simplex K)
        """
        filtration = []
        for ord_simplex in sorted_simplices_list:
            new_filtration = str(len(ord_simplex) - 1) + str(' ') + str(len(ord_simplex) - 1) + str(' ') + "".join(
                str(ord_simplex[i]) + str(' ') for i in range(len(ord_simplex)))
            filtration.append(new_filtration)
        return filtration

    def compute_filtrations(self, dimensions):
        # Compute all the required filtrations
        """
        Input: wanted dimension for the d-ball and for the d-sphere
        Output: list with all the required filtrations, according to the syntax given in "from simplices to filtration"
        """
        all_filtrations = []

        # For the d-ball
        for d in dimensions:
            self.all_names.append("{}-ball".format(d))
            all_filtrations.append(self.from_simplices_to_filtration(
                self.all_simplices_from_top_simplices(self.subsimplices_from_simplex(self.top_simplex_d_ball(d)))))

            # For the d-sphere
            self.all_names.append("{}-sphere".format(d))
            all_filtrations.append(self.from_simplices_to_filtration(
                self.all_simplices_from_top_simplices(self.subsimplices_from_simplex(self.top_simplex_d_ball(d))))[
                                   :-1])  # We remove the (d+1) simplex

        # For the Moebius band
        self.all_names.append("moebius_band")
        all_filtrations.append(self.from_simplices_to_filtration(
            self.all_simplices_from_top_simplices(self.top_simplices_from_figure(self.triangulation_Moebius_band()))))

        # For the torus
        self.all_names.append("torus")
        all_filtrations.append(self.from_simplices_to_filtration(
            self.all_simplices_from_top_simplices(self.top_simplices_from_figure(self.triangulation_torus()))))

        # For the Klein bottle
        self.all_names.append("klein_bottle")
        all_filtrations.append(self.from_simplices_to_filtration(
            self.all_simplices_from_top_simplices(self.top_simplices_from_figure(self.triangulation_Klein_bottle()))))

        # For the projective plane
        self.all_names.append("projective_plane")
        all_filtrations.append(
            self.from_simplices_to_filtration(self.all_simplices_from_top_simplices(self.top_simplices_projective_plane())))

        self.all_filtrations = all_filtrations

    def triangulation_Moebius_band(self):
        """
        Output: 2D-list from the vertices, row by row, corresponding to a triangulation from the Moebius band
        """
        return [['A', 'B', 'C', 'D'], ['D', 'E', 'F', 'A']]

    def triangulation_torus(self):
        """
        Output: 2D-list from the vertices, row by row, corresponding to a triangulation from the torus
        """
        return [['A', 'B', 'C', 'A'], ['D', 'E', 'F', 'D'], ['G', 'H', 'I', 'G'], ['A', 'B', 'C', 'A']]

    def triangulation_Klein_bottle(self):
        """
        Output: 2D-list from the vertices, row by row, corresponding to a triangulation from the Klein bottle
        """
        return [['A', 'B', 'C', 'A'], ['D', 'E', 'F', 'G'], ['G', 'H', 'I', 'D'], ['A', 'B', 'C', 'A']]

    def top_simplices_projective_plane(self):
        """
        Output: set from the top-dimensional simplices (here the 2D simplices from a triangulation from the 2D projective plane), i.e. set of list of vertices
        """
        return set(
            [('A', 'B', 'D'), ('B', 'C', 'D'), ('C', 'D', 'E'), ('C', 'E', 'A'), ('A', 'E', 'B'), ('B', 'E', 'F'),
             ('B', 'C', 'F'), ('A', 'F', 'C'), ('A', 'F', 'D'), ('D', 'E', 'F')])

    def top_simplex_d_ball(self, d):
        """
        Input: the dimension d
        Output: tuple of the (d+1) vertices from the top-dimensional simplex of the d-ball
        """
        assert d < 11  # We could have gone further, one can easily modify it by extending list_max_vertices
        list_vertices = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']
        return tuple(list_vertices[:d + 1])

    # For practical reasons for the drawing of the filtration,
    # we define another triangulation for the Moebius band, the torus, the Klein bottle
    # as the same vertices appear several times

    def triangulation_Moebius_band_for_drawing(self):
        """
        Output: 2D-list from the vertices, row by row, corresponding to a triangulation from the Moebius band
        """
        return [['A1', 'B', 'C', 'D1'], ['D2', 'E', 'F', 'A2']]

    def triangulation_torus_for_drawing(self):
        """
        Output: 2D-list from the vertices, row by row, corresponding to a triangulation from the torus where we 'split' each vertex.
        """
        return [['A1', 'B1', 'C1', 'A2'], ['D1', 'E', 'F', 'D2'], ['G1', 'H', 'I', 'G2'], ['A3', 'B2', 'C2', 'A4']]

    def triangulation_Klein_bottle_for_drawing(self):
        """
        Output: 2D-list from the vertices, row by row, corresponding to a triangulation from the Klein bottle
        """
        return [['A1', 'B1', 'C1', 'A2'], ['D1', 'E', 'F', 'G1'], ['G2', 'H', 'I', 'D2'], ['A3', 'B2', 'C2', 'A4']]

    def compute_all_simplices_for_drawing(self, d):
        """
        Input: wanted dimension for the d-ball and for the d-sphere
        Output: list with all the required simplicial complexes with splitted vertices
        """
        all_filtrations = []

        # For the d-ball
        all_filtrations.append(self.all_simplices_from_top_simplices(self.subsimplices_from_simplex(self.top_simplex_d_ball(d))))

        # For the d-sphere
        all_filtrations.append(self.all_simplices_from_top_simplices(self.subsimplices_from_simplex(self.top_simplex_d_ball(d)))[
                               :-1])  # We remove the (d+1) simplex

        # For the Moebius band
        all_filtrations.append(
            self.all_simplices_from_top_simplices(self.top_simplices_from_figure(self.triangulation_Moebius_band_for_drawing())))

        # For the torus
        all_filtrations.append(
            self.all_simplices_from_top_simplices(self.top_simplices_from_figure(self.triangulation_torus_for_drawing())))

        # For the Klein bottle
        all_filtrations.append(
            self.all_simplices_from_top_simplices(self.top_simplices_from_figure(self.triangulation_Klein_bottle_for_drawing())))

        # For the projective plane
        all_filtrations.append(self.all_simplices_from_top_simplices(self.top_simplices_projective_plane()))

        return all_filtrations

    def filtrations_to_file(self, header):
        for i, name in enumerate(self.all_names):
            filtration = self.all_filtrations[i]
            output_file = open(header + name + ".txt", "w")
            output_file.write("\n".join(filtration))
            output_file.close()

