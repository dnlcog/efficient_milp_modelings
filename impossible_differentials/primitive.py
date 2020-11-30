from gurobipy import *
import pickle
import utilities
from itertools import product as itp
import random
import time


class Primitive:
    """ Gurobi Model of a cryptographic primitive for differential properties. """

    def __init__(self, in_size, out_size):
        # Dictionnary of Sbx modelings used in this primitive.
        self.sbox_modelings = {}

        # Input and output sizes.
        self.in_size = in_size
        self.out_size = out_size

        # Input and output Gurobi variables.
        self.in_var = {}
        self.out_var = {}

        # Miscellaneous objects
        self.misc = {}

        # Gurobi Model
        self.model = Model()

    def add_sbox_modeling(self, file_name, other_name=None):
        """
        Loads the pickle file file_name of an Sbox modeling.
        The pickle file should contain a set of lists representing
        inequalities.
        If ineg is such a list, the input coefficients are first,
        then come the output coefficients and finally comes the constant.
        The inequality is then:
        sum(input[i] * ineg[i]) + sum(output[i] * ineg[i + len(input)])
        + ineg[len(input) + len(output)] >= 0
        """
        with open(file_name, "rb") as f:
            if other_name is None:
                other_name = file_name
            (in_size, out_size, ddt, ineq) = pickle.load(f)
            self.sbox_modelings[other_name] = (
                utilities.ddt_rows(ddt, in_size, out_size),
                utilities.ddt_cols(ddt, in_size, out_size),
                ineq,
            )

    def add_sbox_constr(self, sbox_name, a, b):
        """ 
        Adds constraints for one sbox registered in 
        add_sbox_modelings[sbox_name].
        a is a list of input variables,
        b is a list of output variables,
        """

        n = len(a)
        m = len(b)
        (_, _, ineqs) = self.sbox_modelings[sbox_name]
        # ineqs = self.sbox_modelings[sbox_name]
        for ineg in ineqs:
            assert len(ineg) == n + m + 1
            self.model.addConstr(
                quicksum(ineg[i] * a[i] for i in range(n))
                + quicksum(ineg[i + n] * b[i] for i in range(m))
                + ineg[n + m]
                >= 0
            )

    def add_xor_constr(self, variables, offset=0, mode="binary"):
        """
        If mode = "binary", adds the $2^{n-1}$ constraints modeling
        the XOR constraint x[0] ^ ... ^ x[n-1] = offset
        where x is variables.
        If mode = "integer", models the same XOR constraint with a dummy
        integer variable t with x[0] + ... + x[n-1] = 2 * t + offset.
        """
        x = variables
        n = len(x)

        if mode == "binary" or mode == "both":
            for i in range(1 << n):
                bit_list = utilities.bits(i, n)
                if sum(bit_list) % 2 == (1 - offset):
                    constraint = quicksum(
                        x[j] if bit_list[j] == 0 else 1 - x[j] for j in range(n)
                    )
                    self.model.addConstr(constraint >= 1)
        if mode == "integer" or mode == "both":
            offset = offset % 2

            t = self.model.addVar(
                name="dummy_xor", lb=0, ub=(n // 2) + (n % 2), vtype=GRB.INTEGER
            )
            self.model.addConstr(quicksum(x) == (2 * t) + offset)

    def add_bin_matrix_constr(self, matrix, x, b, mode="binary"):
        """
        Adds constraints given by matrix * x = b
        where x is a list of GRB.BINARY variables
        and b is a constant given as an integer.
        """
        y = utilities.bits(b, len(matrix))

        for i in range(len(matrix)):
            row = matrix[i]
            assert len(row) == len(x)
            variables = [x[j] for j in range(len(x)) if row[j] != 0]
            self.add_xor_constr(variables, offset=y[i], mode=mode)

    def last_input_diff(self):
        """
        Returns the value of the input difference in the last solution.
        """
        x = 0
        for i in range(self.in_size):
            if self.in_var[i].x >= 0.5:
                x ^= 1 << i

        return x

    def last_output_diff(self):
        """
        Returns the value of the output difference in the last solution.
        """
        x = 0
        for i in range(self.out_size):
            if self.out_var[i].x >= 0.5:
                x ^= 1 << i

        return x

    def set_input_diff(self, in_diff):
        """
        Sets the input difference for impossible differential
        search.
        """
        bit_list = utilities.bits(in_diff, self.in_size)

        for i in range(self.in_size):
            key = "in_{}".format(i)
            if key in self.misc:
                self.model.remove(self.misc[key])
            self.misc[key] = self.model.addConstr(self.in_var[i] == bit_list[i])

    def set_output_diff(self, out_diff):
        """
        Sets the output difference for impossible differential
        search.
        """
        bit_list = utilities.bits(out_diff, self.out_size)

        for i in range(self.out_size):
            key = "out_{}".format(i)
            if key in self.misc:
                self.model.remove(self.misc[key])
            self.misc[key] = self.model.addConstr(self.out_var[i] == bit_list[i])

    def set_search_space(self, the_set):
        """
        Sets the search space for the next search of
        impossible differentials.
        """
        self.misc["imp_diff_search_space"] = the_set

    def format_state(self, x):
        """
        Gives a nice representation of a state.
        """
        return "{}".format(x)

    def is_possible(self, x, y):
        """
        Outputs whether the pair (x, y) is a possible transition
        or not.
        """
        self.set_input_diff(x)
        self.set_output_diff(y)
        self.model.optimize()
        status = self.model.status
        output_choices = [
            gurobipy.GRB.OPTIMAL,
            gurobipy.GRB.INFEASIBLE,
        ]
        assert status in output_choices
        return status == gurobipy.GRB.OPTIMAL


class AesLike(Primitive):
    """
    Gurobi Model of an AES-like primitive.
    Ie, when there is only one permutation SBox
    and one linear layer,
    """

    def __init__(self, state_size, nb_rounds, sbox_file):

        Primitive.__init__(self, state_size, state_size)

        self.nb_rounds = nb_rounds
        self.sbox_name = sbox_file

        with open(sbox_file, "rb") as f:
            (in_nibble_size, out_nibble_size, _, _) = pickle.load(f)

        assert in_nibble_size == out_nibble_size
        nibble_size = in_nibble_size

        assert state_size % nibble_size == 0
        nb_nibbles = state_size // nibble_size

        self.add_sbox_modeling(sbox_file)

        self.state_size = state_size
        self.nibble_size = nibble_size
        self.nb_nibbles = nb_nibbles

        in_sbox = {}  # in_sbox  for sbox input
        out_sbox = {}  # out_sbox for sbox output

        for i, j in itp(range(nb_rounds), range(state_size)):
            in_sbox[i, j] = self.model.addVar(
                name="in_sbox_\{%s, %s\}" % (i, j), vtype=gurobipy.GRB.BINARY
            )
        for i, j in itp(range(nb_rounds), range(state_size)):
            out_sbox[i, j] = self.model.addVar(
                name="out_sbox_\{%s, %s\}" % (i, j), vtype=gurobipy.GRB.BINARY
            )

        for i in range(state_size):
            self.in_var[i] = in_sbox[0, i]
            self.out_var[i] = out_sbox[nb_rounds - 1, i]

        for i in range(nb_rounds):
            self.subcell(
                [in_sbox[i, j] for j in range(state_size)],
                [out_sbox[i, j] for j in range(state_size)],
            )

        for i in range(nb_rounds - 1):
            self.linear_layer(
                [out_sbox[i, j] for j in range(state_size)],
                [in_sbox[i + 1, j] for j in range(state_size)],
            )

        self.in_sbox = in_sbox
        self.out_sbox = out_sbox

        random.seed()

    def subcell(self, in_sbox, out_sbox):
        n = self.nb_nibbles
        d = self.nibble_size
        for nibble in range(n):
            a = [in_sbox[(d * nibble) + i] for i in range(d)]
            b = [out_sbox[(d * nibble) + i] for i in range(d)]
            self.add_sbox_constr(self.sbox_name, a, b)

    def linear_layer(self, x_in, x_out):
        """
        Adds linear layer constraints for one round.
        """
        # To be implemented for each primitive
        raise NotImplementedError

    def format_state(self, x):
        """
        Gives a nice representation of a state.
        A state is represented as a list of tuples
        with the cell number and the value of this cell.
        """
        if x == 0:
            return "0"

        string = "["
        bit_list = utilities.bits(x, self.state_size)
        for i in range(self.state_size):
            d = self.nibble_size
            bit = i % d

            if bit == 0:
                val = 0
            val ^= bit_list[i] << bit

            if bit == d - 1:
                if val != 0:
                    string += "(nib: {:2}, val: {:3x}), ".format(i // d, val)

        return string[:-2] + "]"

    def get_state_in_sbox(self, r):
        """
        Gets the state at the input of round r.
        """
        assert r < self.nb_rounds
        out = 0
        for i in range(self.state_size):
            if self.in_sbox[r, i].x >= 0.5:
                out ^= 1 << i
        return out

    def get_state_out_sbox(self, r):
        """
        Gets the state at the output of round r.
        """
        assert r < self.nb_rounds
        out = 0
        for i in range(self.state_size):
            if self.out_sbox[r, i].x >= 0.5:
                out ^= 1 << i
        return out

    def equimip_search(self, the_dict, aux_in, aux_out, message=""):
        """
        More general version of the differential possibility equivalence technique
        of Sasaki and Todo EC17.
        the_dict: map (python dict) from an input difference to try to
            a set of output differences to try with it.
        aux_in: auxiliary input model of the same class with a smaller
            number of rounds.
        aux_out: same for output.
        """

        r_in = aux_in.nb_rounds - 1
        r_out = aux_out.nb_rounds - 1
        assert r_in + r_out < self.nb_rounds
        assert r_in >= 0
        assert r_out >= 0

        # Useful for printing progress. START
        def remaining(the_dict):
            return sum([len(the_dict[x]) for x in iter(the_dict)])

        def spaces(x):
            return "".join([" " for i in range(x)])

        out = []
        length = remaining(the_dict)
        nb_done = 0
        nb_milp = 0
        nb_milp_x = 0
        nb_milp_y = 0
        discarded = 0

        absolute_time = time.time()

        print(
            "| {}Message |".format(spaces(len(message) - 7))
            + " {}Time |".format(spaces(17 - 4))
            + "  Results / {:10} |".format(length)
            + " {}MIP queries |".format(spaces(20 - 11))
            + "  x queries |"
            + "  y queries |"
            + " Dis. rate |"
            + " Found |"
        )

        def outer_printer():
            seconds = int(time.time() - absolute_time)
            minutes = seconds // 60
            hours = minutes // 60
            str_time = "{:4}h, {:2}min, {:2}s".format(hours, minutes % 60, seconds % 60)

            return (
                "| {} | {} |  {:5.1f} % = {:10} |".format(
                    message, str_time, (100.0 * nb_done) / length, nb_done,
                )
                + " {:10} = {:5.2f} % |".format(nb_milp, (100.0 * nb_milp) / length,)
                + " {:10} |".format(nb_milp_x,)
                + " {:10} |".format(nb_milp_y,)
                + "   {:5.1f} % |".format(
                    (100.0 * discarded) / (nb_milp_x + nb_milp_y)
                    if nb_milp_x + nb_milp_y != 0
                    else 0,
                )
                + " {:5} |".format(len(out))
            )

        # END

        # While there are difference pairs to try...
        while len(the_dict) >= 1:
            # x is the input difference we are going to try.
            x = list(the_dict)[0]
            while len(the_dict[x]) >= 1:
                outer_print = outer_printer()
                print(outer_print, end="\n")

                # y is the input difference we are going to try.
                y = the_dict[x].pop()

                # Printing this message while the solver is running
                # on the main model self (in the function self.is_possible)
                # This computation can last for a few hours.
                print(
                    "MIP query on input {} and output {}".format(
                        self.format_state(x), self.format_state(y),
                    ),
                    end="\r",
                )
                possible = self.is_possible(x, y)
                nb_milp += 1

                # If we have found an impossible differential, add it to the output.
                if not possible:
                    out.append((x, y))
                # Else use the differential possibility equivalence technique.
                else:
                    # Get the middle values in the computed path.
                    x_mid = self.get_state_out_sbox(r_in)
                    y_mid = self.get_state_in_sbox(self.nb_rounds - r_out - 1)

                    to_discard = []

                    # Printer related stuff. START

                    visited = 0
                    rem = remaining(the_dict)

                    class inner_printer:
                        def __init__(self):
                            self.timer = time.time()

                        def go(self):
                            if (time.time() - self.timer) >= 0.5:
                                self.timer = time.time()
                                print(
                                    "Discarding progress "
                                    + "{:.1f} %   rate {:.1f} %".format(
                                        (100.0 * visited) / rem,
                                        (100.0 * len(to_discard)) / visited,
                                    )
                                    + spaces(30),
                                    end="\r",
                                )

                    ip = inner_printer()

                    # END

                    # For each possible input difference...
                    for x_start in the_dict.keys():
                        # We first try to compute the beginning of the path
                        # between x_start and x_mid (if x_start is not the initial x).
                        try_y = x == x_start
                        if not try_y:
                            nb_milp_x += 1
                            try_y = aux_in.is_possible(x_start, x_mid)

                        # If a path from x_start to x_mid is found...
                        if try_y:
                            # For each output difference y to try with input x_start...
                            for y in the_dict[x_start]:
                                nb_milp_y += 1
                                visited += 1
                                # We check whether there is a path between y_mid and y.
                                if aux_out.is_possible(y_mid, y):
                                    # If it is the case, we will discard the pair (x_start, y)
                                    # from input/output pairs to try.
                                    to_discard.append((x_start, y))
                                    discarded += 1
                                ip.go()
                        else:
                            visited += len(the_dict[x_start])
                            ip.go()

                    for (x, y) in to_discard:
                        the_dict[x].remove(y)

                nb_done = length - remaining(the_dict)

            # Check and clean the set of input/output pairs to try.
            assert len(the_dict[x]) == 0
            keys = list(the_dict)
            for key in keys:
                if len(the_dict[key]) == 0:
                    del the_dict[key]

        nb_done = length - remaining(the_dict)
        outer_print = outer_printer()
        print(outer_print, end="\n")

        return out

    def minimize_active_sboxes(self):
        """
        Computes the minimum number of active SBoxes.
        """
        local_constraints = []

        # variables in y model whether a CELL is active or not.
        y = dict()
        for r, i in itp(range(self.nb_rounds), range(self.nb_nibbles)):
            y[r, i] = self.model.addVar(
                name="active_({}, {})".format(r, i), vtype=GRB.BINARY, obj=1.0,
            )
            bits = [
                self.in_sbox[r, (self.nibble_size * i) + j]
                for j in range(self.nibble_size)
            ]
            constr = self.model.addGenConstrOr(y[r, i], bits)
            local_constraints.append(constr)

        # We fix at least one active input cell.
        constr = self.model.addConstr(
            quicksum(y[0, i] for i in range(self.nb_nibbles)) >= 1
        )
        local_constraints.append(constr)

        self.model.optimize()

        output = []
        if self.model.status == gurobipy.GRB.OPTIMAL:
            for r, i in itp(range(self.nb_rounds), range(self.nb_nibbles)):
                if y[r, i].x >= 0.5:
                    output.append((r, i))
            assert len(output) == int(self.model.objVal)

        else:
            print("Optimization ended with status {}".format(self.model.status))
            exit(1)

        # We clean the model from local constraints and variables.
        for constr in local_constraints:
            self.model.remove(constr)
        for r, i in itp(range(self.nb_rounds), range(self.nb_nibbles)):
            self.model.remove(y[r, i])

        return output
