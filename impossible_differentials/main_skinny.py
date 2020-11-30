from skinny import *
import argparse
import itertools

if __name__ == "__main__":

    nb_rounds = 13

    parser = argparse.ArgumentParser(
        description="Launches search for impossible differentials on 13 rounds of Skinny "
        + "with chosen active input cell."
    )
    parser.add_argument(
        "cell", type=int, help="Active input cell.", choices=[i for i in range(16)],
    )
    args = parser.parse_args()

    cell = args.cell

    # Main model
    mid = Skinny(nb_rounds, "skinny_sbox.pkl")
    mid.model.setParam("LogToConsole", 0)

    # Auxiliary input model
    aux_in = Skinny(2, "skinny_sbox.pkl")
    aux_in.model.setParam("LogToConsole", 0)

    # Auxiliary output model
    aux_out = Skinny(2, "skinny_sbox.pkl")
    aux_out.model.setParam("LogToConsole", 0)

    for out_cell in range(16):
        the_dict = dict()
        for valx in range(1, 1 << 8):
            x = valx << (8 * cell)
            the_dict[x] = set()
            for valy in range(1, 1 << 8):
                y = valy << (8 * out_cell)
                the_dict[x].add(y)

        message = "Skinny {}r in {} out {}.".format(nb_rounds, cell, out_cell)
        res = mid.equimip_search(the_dict, aux_in, aux_out, message=message)
