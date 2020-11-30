import gurobipy
from primitive import AesLike
from itertools import product as itp
import utilities

# Compute the ShiftRows transformation.
shift_rows = [0] * 16
for nib in range(16):
    row = nib % 4
    new_col = nib // 4
    old_col = (new_col + row) % 4
    shift_rows[nib] = (4 * old_col) + row

# A = (M|I) where M is the MixColumns matrix. Little Endian.
mixcolA_original = [
    0x101018180,
    0x202028381,
    0x404040602,
    0x808088C84,
    0x1010109888,
    0x2020203010,
    0x4040406020,
    0x808080C040,
    0x10001818001,
    0x20002838102,
    0x40004060204,
    0x800088C8408,
    0x100010988810,
    0x200020301020,
    0x400040602040,
    0x800080C04080,
    0x1000081800101,
    0x2000083810202,
    0x4000006020404,
    0x800008C840808,
    0x10000098881010,
    0x20000030102020,
    0x40000060204040,
    0x800000C0408080,
    0x100000080010181,
    0x200000081020283,
    0x400000002040406,
    0x80000008408088C,
    0x1000000088101098,
    0x2000000010202030,
    0x4000000020404060,
    0x80000000408080C0,
]

# A = P x (M|I) x Q. Matrix that needs less inequalities.
mixcolA = [
    0x1010101000080,
    0x202020301,
    0x404040602,
    0x8000808840C0000,
    0x1010101808,
    0x2020203010,
    0x4040406020,
    0x808080C040,
    0x101000100800100,
    0x20002030102,
    0x40004060204,
    0x80808000000840C,
    0x100010180810,
    0x200020301020,
    0x400040602040,
    0x800080C04080,
    0x101010000008001,
    0x2000003010202,
    0x4000006020404,
    0x808080C000084,
    0x10000018081010,
    0x20000030102020,
    0x40000060204040,
    0x800000C0408080,
    0x100010180010000,
    0x200000001020203,
    0x400000002040406,
    0x808000800840C00,
    0x1000000008101018,
    0x2000000010202030,
    0x4000000020404060,
    0x80000000408080C0,
]

# A = (I|I). Useful for testing purposes.
mixcolA_identity = [(1 ^ (1 << 32)) << i for i in range(32)]

# Those functions apply the affine change in the SBox when we choose
# the affine equivalent SBox.
def qmat(x):
    y = x ^ ((x >> 7) & 1) ^ (((x >> 7) & 1) << 3)
    return y


def qmat_on_state(x):
    out = 0
    for i in range(16):
        out ^= qmat((x >> (8 * i)) & 0xFF) << (8 * i)
    return out


class Aes(AesLike):
    """ Gurobi Model for AES differential trails. """

    def __init__(self, nb_rounds, sbox_file, mixcol="equiv"):

        # Different choices of MixColumns models.
        if mixcol == "equiv":
            mixcol_hex = mixcolA
        elif mixcol == "origin":
            mixcol_hex = mixcolA_original
        else:
            mixcol_hex = mixcolA_identity

        mixcol_matrix = []
        for row in mixcol_hex:
            bit_list = utilities.bits(row, 64)
            mixcol_matrix += [bit_list]

        self.mat = mixcol_matrix
        self.mixcol = mixcol

        AesLike.__init__(self, 128, nb_rounds, sbox_file)

    def linear_layer(self, x_in, x_out):
        # Reversing the bytes because of the AES bytes being
        # ordered in big endian.
        x_in2 = []
        for nib in range(16):
            x_in2 += [x_in[(8 * (15 - nib)) + i] for i in range(8)]
        x_in2 = x_in

        x_out2 = []
        for nib in range(16):
            x_out2 += [x_out[(8 * (15 - nib)) + i] for i in range(8)]
        x_out2 = x_out

        x_in = [
            x_in[(8 * shift_rows[i // 8]) + (i % 8)] for i in range(128)
        ]  # variables after the shift_rows, at the input of mixcol

        for col in range(4):
            bit_list = [x_in[(32 * col) + i] for i in range(32)] + [
                x_out[(32 * col) + i] for i in range(32)
            ]

            self.add_bin_matrix_constr(self.mat, bit_list, 0, mode="binary")

    def set_output_diff(self, out_diff):
        """
        Sets the output difference for impossible differential
        search. Replaces the original function if this model
        uses the equivalent linear mixcolumns. Indeed, this
        mixcolumns matrix has a modified input and then a
        modified output for the Sbox. This function makes it
        transparent.
        """
        if self.mixcol == "equiv":
            x = qmat_on_state(out_diff)
        else:
            x = out_diff
        AesLike.set_output_diff(self, x)

    def last_output_diff(self):
        """
        Returns the value of the output difference in the last solution.
        Replacement for the same reason as set_output_diff.
        """
        out_diff = AesLike.last_output_diff(self)

        if self.mixcol == "equiv":
            return qmat_on_state(out_diff)
        else:
            return out_diff

    def get_state_out_sbox(self, r):
        """
        Gets the state at the output of round r.
        Replacement for the same reason as set_output_diff.
        """
        x = AesLike.get_state_out_sbox(self, r)
        if self.mixcol == "equiv":
            x = qmat_on_state(x)

        return x

    def set_active_input_cell(self, cell):
        """
        Adds constraints that speed up the computation when one input byte is active.
        Those new constraints just ecplicit which variables are 0 for sure in and out
        the first and second SBox layers.
        """
        assert 0 <= cell and cell < 16
        # Checks that no active input cell has been set before.
        for i in range(16):
            key = "in_cell_{}".format(i)
            assert key not in self.misc

        key = "in_cell_{}".format(cell)
        self.misc[key] = set()

        # For each state bit...
        for i in range(128):
            # Fix 0 values in and out the first SBox layer.
            if i // 8 != cell:
                self.misc[key].add(self.model.addConstr(self.in_sbox[0, i] == 0))
                self.misc[key].add(self.model.addConstr(self.out_sbox[0, i] == 0))

            # Fix 0 values in and out the second SBox layer.
            if self.nb_rounds >= 2 and shift_rows[i // 8] // 4 != cell // 4:
                self.misc[key].add(self.model.addConstr(self.in_sbox[1, i] == 0))
                self.misc[key].add(self.model.addConstr(self.out_sbox[1, i] == 0))

    def unset_active_input_cell(self):
        """
        Removes active input bytes constraints.
        """
        for cell in range(16):
            key = "in_cell_{}".format(cell)
            if key in self.misc:
                for constr in self.misc[key]:
                    self.model.remove(constr)
                del self.misc[key]

    def unset_active_output_cell(self):
        """
        Same as above but for the output.
        """
        for cell in range(16):
            key = "out_cell_{}".format(cell)
            if key in self.misc:
                for constr in self.misc[key]:
                    self.model.remove(constr)
                del self.misc[key]

    def set_active_output_cell(self, cell):
        """
        Same as above but for the output.
        """
        assert 0 <= cell and cell < 16
        for i in range(16):
            key = "out_cell_{}".format(i)
            assert key not in self.misc

        key = "out_cell_{}".format(cell)
        self.misc[key] = set()

        for i in range(128):
            if i // 8 != cell:
                self.misc[key].add(
                    self.model.addConstr(self.in_sbox[self.nb_rounds - 1, i] == 0)
                )
                self.misc[key].add(
                    self.model.addConstr(self.out_sbox[self.nb_rounds - 1, i] == 0)
                )

            if self.nb_rounds >= 2 and (i // 8) // 4 != shift_rows[cell] // 4:
                self.misc[key].add(
                    self.model.addConstr(self.in_sbox[self.nb_rounds - 2, i] == 0)
                )
                self.misc[key].add(
                    self.model.addConstr(self.out_sbox[self.nb_rounds - 2, i] == 0)
                )


# Some test vectors from the FIPS197.
original_tv = [
    (0x84FB386F1AE1AC977941DD70832DD769, 0x9F487F794F955F662AFC86ABD7F1AB29),
    (0x1F770C64F0B579DEAAAC432C3D37CF0E, 0xB7A53ECBBF9D75A0C40EFC79B674CC11),
    (0x684AF5BC0ACCE85564BB0878242ED2ED, 0x7A1E98BDACB6D1141A6944DD06EB2D3E),
    (0x9316DD47C2FA92834390A1DE43E43F23, 0xAAA755B34CFFE57CEF6F98E1F01C13E6),
]

shift_rows_tv = [
    (0x63CAB7040953D051CD60E0E7BA70E18C, 0x6353E08C0960E104CD70B751BACAD0E7),
    (0x84FB386F1AE1AC97DF5CFD237C49946B, 0x84E1FD6B1A5C946FDF4938977CFBAC23),
]

equiv_tv = [
    (0x00102030405060708090A0B0C0D0E0F0, 0x63CAB7040953D051CD60E0E7BA70E18C),
    (0x89D810E8855ACE682D1843D8CB128FE4, 0xA761CA9B97BE8B45D8AD1A611FC97369),
    (0x4915598F55E5D7A0DACA94FA1F0A63F7, 0x3B59CB73FCD90EE05774222DC067FB68),
    (0xFA636A2825B339C940668A3157244D17, 0x2DFB02343F6D12DD09337EC75B36E3F0),
    (0x247240236966B3FA6ED2753288425B6C, 0x36400926F9336D2D9FB59D23C42C3950),
]

# For stupid endianness reasons...
def reverse(x):
    y = 0
    for i in range(16):
        y ^= ((x >> (8 * i)) & 0xFF) << (8 * (15 - i))
    return y


original_tv = [(reverse(x), reverse(y)) for (x, y) in original_tv]
shift_rows_tv = [(reverse(x), reverse(y)) for (x, y) in shift_rows_tv]
equiv_tv = [(reverse(x), reverse(y)) for (x, y) in equiv_tv]


def test_shift_rows():
    mid = Aes(2, "identity_sbox_8.pkl", mixcol="identity")
    mid.model.setParam("LogToConsole", 0)
    for (x, y) in shift_rows_tv:
        mid.set_input_diff(x)
        mid.set_output_diff(y)
        mid.model.optimize()
        assert mid.model.status == gurobipy.GRB.OPTIMAL
    print("SR ok")


def test_original_lin_layer():
    mid = Aes(2, "identity_sbox_8.pkl", mixcol="origin")
    mid.model.setParam("LogToConsole", 0)
    for (x, y) in original_tv:
        mid.set_input_diff(x)
        mid.set_output_diff(y)
        mid.model.optimize()
        assert mid.model.status == gurobipy.GRB.OPTIMAL

        mid.set_input_diff(x ^ 1)
        mid.model.optimize()
        assert mid.model.status == gurobipy.GRB.INFEASIBLE
    print("Original Linear Layer ok")


def test_equiv_lin_layer():
    mid = Aes(2, "identity_sbox_8.pkl", mixcol="equiv")
    mid.model.setParam("LogToConsole", 0)
    for (x, y) in original_tv:
        xx = 0
        for i in range(16):
            xx ^= qmat((x >> (8 * i)) & 0xFF) << (8 * i)
        yy = 0
        for i in range(16):
            yy ^= qmat((y >> (8 * i)) & 0xFF) << (8 * i)
        mid.set_input_diff(xx)
        mid.set_output_diff(yy)
        mid.model.optimize()
        assert mid.model.status == gurobipy.GRB.OPTIMAL

        mid.set_input_diff(xx ^ 1)
        mid.model.optimize()
        assert mid.model.status == gurobipy.GRB.INFEASIBLE
    print("Original Equiv Linear Layer ok")


def test_equiv_sbox(sbox_file):
    mid = Aes(1, sbox_file, mixcol="equiv")
    mid.model.setParam("LogToConsole", 0)
    nb_vec = len(equiv_tv) // 2
    test_vectors = [
        (x0 ^ x1, y0 ^ y1)
        for ((x0, y0), (x1, y1)) in itp(equiv_tv, equiv_tv)
        if x0 != x1
    ]
    xx = 0
    for i in range(16):
        xx ^= i << (8 * i)
    for (x, y) in test_vectors:
        mid.set_input_diff(x)
        mid.set_output_diff(y)
        mid.model.optimize()
        assert mid.model.status == gurobipy.GRB.OPTIMAL
        mid.set_input_diff(x ^ xx)
        mid.model.optimize()
        assert mid.model.status == gurobipy.GRB.INFEASIBLE
    print("Equiv Sbox ok")


if __name__ == "__main__":
    test_shift_rows()
    test_original_lin_layer()
    test_equiv_lin_layer()
    test_equiv_sbox("greedy_sbox_ineg_aes_equiv.pkl")
