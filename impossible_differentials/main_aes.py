from aes import *
import argparse
import itertools

if __name__ == "__main__":

    nb_rounds = 5

    parser = argparse.ArgumentParser(
        description="Launches search for impossible differentials in 5 rounds of AES "
        + "with chosen active input and output cells."
    )
    parser.add_argument(
        "in_cell", type=int, help="Nibble of input.", choices=[i for i in range(16)],
    )
    parser.add_argument(
        "out_cell", type=int, help="Nibble of output.", choices=[i for i in range(16)],
    )
    args = parser.parse_args()

    in_cell = args.in_cell
    out_cell = args.out_cell

    # Main model.
    mid = Aes(nb_rounds, "aes_equiv_sbox.pkl")
    mid.model.setParam("LogToConsole", 0)
    mid.set_active_input_cell(in_cell)
    mid.set_active_output_cell(out_cell)

    # Auxiliary input model.
    aux_in = Aes(2, "aes_equiv_sbox.pkl")
    aux_in.model.setParam("LogToConsole", 0)
    aux_in.set_active_input_cell(in_cell)

    # Auxiliary output model.
    aux_out = Aes(2, "aes_equiv_sbox.pkl")
    aux_out.model.setParam("LogToConsole", 0)
    aux_out.set_active_output_cell(out_cell)

    the_dict = dict()
    for valx in range(1, 1 << 8):
        x = valx << (8 * in_cell)
        the_dict[x] = set()
        for valy in range(1, 1 << 8):
            y = valy << (8 * out_cell)
            the_dict[x].add(y)

    message = "Aes 5r in {} out {}.".format(in_cell, out_cell)
    mid.equimip_search(the_dict, aux_in, aux_out, message=message)
