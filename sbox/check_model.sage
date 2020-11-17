import argparse
import pickle
from sage.crypto.sboxes import SBox
from sage.crypto.sboxes import sboxes
from operator import xor


def point_kept(in_size, out_size, inequality, a, b):
    """
    Checks whether inequality is satisfied by point (a, b)
    """
    result = inequality[in_size + out_size]
    for i in range(in_size):
        result += inequality[i] * ((a >> i) & 1)
    for i in range(out_size):
        result += inequality[i + in_size] * ((b >> i) & 1)
    return result >= 0


def check_ddt_ineq(in_size, out_size, ddt, ineqs):
    """
    Checks the sbox model given by ineqs against the ddt.
    """
    ddt_ones = [(a, b) for (a, b) in ddt if ddt[a, b] != 0]
    ddt_zeros = set([(a, b) for (a, b) in ddt if ddt[a, b] == 0])
    length = len(ineqs)
    i = 0
    for ineq in ineqs:
        print("  {} / {}".format(i, length), end="\r")
        for (a, b) in ddt_ones:
            assert point_kept(in_size, out_size, ineq, a, b)
        elim = []
        for (a, b) in ddt_zeros:
            if not point_kept(in_size, out_size, ineq, a, b):
                elim.append((a, b))
        for (a, b) in elim:
            ddt_zeros.remove((a, b))
        i += 1
    assert len(ddt_zeros) == 0


def aes_equiv():
    """
    Generates the affine equivalent Sbox for the AES
    """

    def qmat(x):
        y = xor(x, xor((x >> 7) & 1, ((x >> 7) & 1) << 3))
        return y

    aes = sboxes["AES"]
    value_table = [qmat(aes[x]) for x in range(256)]
    return SBox(value_table)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Checks wether the sbox model in the pickle file is correct."
    )
    parser.add_argument(
        "sbox",
        type=str,
        choices=list(sboxes.keys()) + ["AES_equiv"],
        help="Sbox that should be modeled.",
    )
    parser.add_argument(
        "ineq_file",
        type=str,
        help="Pickle file containing model" + "(in_size, out_size, ddt, ineq_set).",
    )
    args = parser.parse_args()

    if args.sbox == "AES_equiv":
        sbox = aes_equiv()
    else:
        sbox = sboxes[args.sbox]

    i = sbox.input_size()
    o = sbox.output_size()
    ddt_sage = sbox.difference_distribution_table()

    with open(args.ineq_file, "rb") as f:
        data = pickle.load(f)

        for x in list(data):
            print(type(x))

        (in_size, out_size, ddt, ineq_set) = data

    # Checks the first three elements of the tuple.
    assert in_size == i
    assert out_size == o
    for a in range(1 << i):
        for b in range(1 << o):
            assert ddt[a, b] == ddt_sage[a, b]

    check_ddt_ineq(in_size, out_size, ddt, ineq_set)

    print(
        "This model of {} with {} inequalities is correct.".format(
            args.sbox, len(ineq_set)
        )
    )
