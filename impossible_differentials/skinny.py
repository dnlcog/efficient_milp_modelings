import gurobipy
from primitive import AesLike
from itertools import product as itp
from itertools import starmap as itsm
import utilities

shift_rows = [0, 1, 2, 3, 7, 4, 5, 6, 10, 11, 8, 9, 13, 14, 15, 12]

mixcol_equiv = [
    [0, 0, 0, 1, 1, 0, 0, 1],
    [1, 0, 0, 0, 0, 1, 0, 0],
    [0, 1, 1, 0, 0, 0, 1, 0],
    [1, 0, 1, 0, 0, 0, 0, 1],
]

mixcol_origin = [
    [1, 0, 1, 1, 1, 0, 0, 0],
    [1, 0, 0, 0, 0, 1, 0, 0],
    [0, 1, 1, 0, 0, 0, 1, 0],
    [1, 0, 1, 0, 0, 0, 0, 1],
]


class Skinny(AesLike):
    """ Gurobi Model for Skinny-128 differential trails. """

    def __init__(self, nb_rounds, sbox_file, mixcol="equiv"):

        # Different choices of MixColumns models.
        if mixcol == "equiv":
            self.mixcol = mixcol_equiv
        else:
            self.mixcol = mixcol_origin
        AesLike.__init__(self, 128, nb_rounds, sbox_file)

    def linear_layer(self, x_in, x_out):
        x_in = [
            x_in[(8 * shift_rows[i // 8]) + (i % 8)] for i in range(128)
        ]  # variables after the shift_rows, at the input of mixcol

        for col, bit in itp(range(4), range(8)):
            bit_list = [x_in[(32 * row) + (8 * col) + bit] for row in range(4)] + [
                x_out[(32 * row) + (8 * col) + bit] for row in range(4)
            ]

            self.add_bin_matrix_constr(
                self.mixcol, bit_list, 0, mode="binary",
            )


def lin_layer(x):
    """ Skinny linear layer implementation for testing purposes. """
    # Collect in cells
    nib = {}
    for i in range(16):
        nib[i] = (x >> (8 * i)) & 0xFF

    # Shift rows
    nib2 = {}
    for row in range(4):
        for j in range(4):
            nib2[(4 * row) + ((j + row) % 4)] = nib[(4 * row) + j]

    # Collect in rows
    row = {}
    for i in range(4):
        row[i] = 0
        for j in range(4):
            row[i] ^= nib2[(4 * i) + j] << (8 * j)

    # Mixcolumns
    row[1] ^= row[2]
    row[2] ^= row[0]
    row[3] ^= row[2]

    row2 = {}
    for i in range(4):
        row2[i] = row[(i - 1) % 4]

    # Collect in int
    out = 0
    for i in range(4):
        out ^= row2[i] << (32 * i)

    return out


def test_linear_layer():
    """
    Testing the modeling of the linear layer with identity Sbox.
    """
    n = 10
    mid = Skinny(n, "identity_sbox_8.pkl")
    mid.model.setParam("LogToConsole", 0)

    the_set = set()
    for i in range(16):
        x = 1 << (8 * i)
        y = x
        for j in range(n - 1):
            y = lin_layer(y)
        the_set.add((x, y))

    res = mid.search_impossible_diff(the_set, message="Linear layer test.")
    assert len(res) == 0

    the_set = set()
    for i in range(16):
        x = 1 << (8 * i)
        y = x
        for j in range(n - 1):
            y = lin_layer(y)
        the_set.add((x, y ^ 1))

    res = mid.search_impossible_diff(the_set, message="Linear layer test.")
    assert len(res) == 16

    print("Linear layer test OK.")


def test_paper_single_impossible_diff():
    """
    Testing the model with arbitrary Sbox against
    the truncated impossible differential trail given in the
    eprint paper:
    The SKINNY Family of Block Ciphers and its Low-Latency Variant MANTIS.
    from the authors of Skinny.
    """
    n = 11
    mid = Skinny(n, "arbitrary_sbox_8_8.pkl")
    mid.model.setParam("LogToConsole", 0)

    the_set = set()
    x = 1 << (8 * 12)
    y = 1 << (8 * 7)
    the_set.add((x, y))

    res = mid.search_impossible_diff(
        the_set, message="Paper single impossible differential."
    )
    assert len(res) == 1

    print("Single impossible differential OK.")


def test_paper_all_impossible_diff():
    """
    Same as above with all the impossible differentials.
    """
    n = 12
    mid = Skinny(n, "arbitrary_sbox_8_8.pkl")
    mid.model.setParam("LogToConsole", 0)

    the_set = set()
    for i in range(16):
        x = 1 << (8 * i)
        for j in range(16):
            y = 1 << (8 * j)
            the_set.add((x, y))

    res = mid.search_impossible_diff(the_set, message="Paper impossible differentials.")

    assert len(res) == 12

    print("All impossible differentials test OK.")


if __name__ == "__main__":
    """
    This section aims at testing this MIP model of Skinny
    for impossible differential search.
    """
    test_linear_layer()
    test_paper_single_impossible_diff()
    test_paper_all_impossible_diff()
