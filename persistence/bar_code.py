import time
class BoundaryMatrix:
    def __init__(self, filtration):
        """compute the boundaries"""
        self.filtration = filtration
        # We define low(i) = the maximal non-zero element for the i-th simplex

        # dic_low[k] = list containing all the simplices i whose low(i) is in the k-th row, i.e. corresponds to the k-th simplex
        # dic_boundaries[i] = list containing all the non-zero boundary simplices from the i-th simplex
        # dic_coord[(i, j)] = 1 if Bi,j with B boundary matrix
        dic_low, dic_boundaries, dic_coord = {}, {}, {}
        for i in range(len(filtration.simplex_list)):
            dic_boundaries[i] = set()  # initialization
            dic_low[i] = set()

        for simplex in filtration.simplex_list:

            simplex_vert = tuple(sorted(simplex.vert))
            simplex_coord = filtration.dic_simplex[simplex_vert]  # identifier of the simplex
            boundaries_id = self.compute_boundary(simplex)
            if len(boundaries_id) > 0:
                dic_boundaries[simplex_coord] = boundaries_id
                low_i = max(boundaries_id)
                dic_low[low_i].add(simplex_coord)


        self.dic_boundaries = dic_boundaries
        self.dic_coord = dic_coord
        self.dic_low = dic_low
        self.gaussian_elim_done = False


    def compute_boundary(self, simplex):
        """return the identifiers of the boundaries of the simplex"""
        boundaries_id = set()
        for i in range(simplex.dim + 1):
            new_vert = list(simplex.vert)
            del new_vert[i]
            if len(new_vert) > 0:
                boundaries_id.add(self.filtration.dic_simplex[tuple(sorted(new_vert))])  # new_vert must be sorted because we compare tuple (it's actually already sorted but we sort it for code clarity)
        return boundaries_id

    def sum_boundaries(self, col1, col2, boundaries1, boundaries2):
        """sum the list of boundaries of 2 simplex in Z/2Z
        Complexity O(max(len(boundaries2), len(boundaries2)))"""
        res = []
        for b in boundaries1:
            if not self.dic_coord.get((col2, b)):
                    res.append(b)
        for b in boundaries2:
            if not self.dic_coord.get((col1, b)):
                res.append(b)
        return res

    def gaussian_elim(self):
        print("beginning gaussian elimination")
        t = time.time()
        for i in reversed(range(0, len(self.dic_boundaries.keys()))):
            if i in self.dic_low:
                columns = sorted(self.dic_low[i])
                if columns:  # if there are columns to kill

                    pivot_column = columns[0]
                    for column in columns[1:]:
                        # sum in Z/2Z
                        self.dic_boundaries[column].symmetric_difference_update(self.dic_boundaries[pivot_column])
                        #update dic_low
                        if self.dic_boundaries[column]:
                            # update the pivot
                            new_pivot = max(self.dic_boundaries[column])
                            self.dic_low[new_pivot].add(column)

        print("gaussian elimination: Done !")
        print("took {} seconds".format(time.time() - t))
        self.gaussian_elim_done = True

    def bar_codes(self, threshold=0.05):
        """compute the bar_code of a filtration (should be computed after gaussian_elim)"""
        print("beginning bar codes computing")
        t = time.time()
        if not self.gaussian_elim_done:
            self.gaussian_elim()
        bar_codes = {}
        bar_codes_list = []
        for i, s in enumerate(self.filtration.simplex_list):
            if len(self.dic_boundaries[i]) == 0:  # create a cycle
                bar_codes[i] = None
            else:  # kill a cycle (create a boundary)
                k = max(self.dic_boundaries[i])  # pivot
                if k in bar_codes.keys():
                    bar_codes[k] = i

        for s1_id in bar_codes:
            dim = self.filtration.simplex_list[s1_id].dim
            s2_id = bar_codes[s1_id]
            s1_val = self.filtration.simplex_list[s1_id].val
            if s2_id:  # the cycle as been killed
                s2_val = self.filtration.simplex_list[s2_id].val
                if s2_val - s1_val > threshold:
                    bar_codes_list.append((dim, s1_val, s2_val))
                    print((dim, s1_val, s2_val))
            else:  # the cycle has never been killed
                s2_val = "inf"
                bar_codes_list.append((dim, s1_val, s2_val))
                print((dim, s1_val, s2_val))

        print("bar codes computing: Done !")
        print("took {} seconds".format(time.time() - t))
        print("for {} simplex".format(len(self.filtration.simplex_list)))

        return bar_codes_list, time.time() - t, len(self.filtration.simplex_list)



