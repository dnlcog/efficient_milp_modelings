import pickle
import sys


def generate_arbitrary_sbox(in_size, out_size):
    """
    Generates an arbitrary sbox pickle file as explained in
    Sasaki Todo EC17.
    """
    ineg_set = set()
    for output in range(out_size):
        ineg = ([1] * in_size) + ([0] * (out_size + 1))
        ineg[output + in_size] = -1
        ineg_set.add(tuple(ineg))
    for inp in range(in_size):
        ineg = ([0] * in_size) + ([1] * out_size) + [0]
        ineg[inp] = -1
        ineg_set.add(tuple(ineg))

    ddt = {}
    for a in range(1 << in_size):
        for b in range(1 << out_size):
            ddt[a, b] = 1 if a != b else 0

    with open("arbitrary_sbox_{}_{}.pkl".format(in_size, out_size), "wb") as f:
        pickle.dump((in_size, out_size, ddt, ineg_set), f, 3)


if __name__ == "__main__":
    try:
        input_size = int(sys.argv[1])
    except IndexError:
        raise SystemExit(
            "Usage: {} ".format(sys.argv[0]) + "<input size> <output size if different>"
        )

    try:
        output_size = int(sys.argv[2])
    except IndexError:
        output_size = input_size

    generate_arbitrary_sbox(input_size, output_size)
