import cplex
import numpy as np
import os
import sys
from numpy import exp, log
from scipy.optimize import broyden1, fsolve
from typing import List, Optional

EPS = 1e-10
MIN_EXP = 1e-320


class Instance:

    def __init__(self, path: str):
        # initalizing class variables
        self.m = 0  # number of weapons
        self.n = 0  # number of targets
        self.W = []  # weapons
        self.T = []  # targets
        self.w = []  # target weights
        self.p = {}  # probabilities of destroying a target by a wepaon
        self.mu = []  # availabilities of weapons

        # reading data from text file
        f = open(path, "r")
        lines = f.readlines()
        f.close()

        # reading first line
        self.m, self.n, mu = [int(v) for v in lines[0].split()]
        # asssuming all weapons have the same availability:
        self.mu = [mu] * self.m

        # reading target weights
        for line in lines[1:self.n + 1]:
            self.w.append(int(line))

        # reading probabilities
        for line in lines[self.n + 1:]:
            splitted = line.split()
            i, j, p = int(splitted[0]), int(splitted[1]), float(splitted[2])
            self.p[i, j] = p

        # creating auxiliary sets
        self.W = list(range(self.m))
        self.T = list(range(self.n))

    def R(self, i: int):
        """
        This function provides a straightforward "syntax sugar" to simplify the code.
        Calling R(i) will provide an iterable [1, ..., mu[i]]
        """
        return range(1, self.mu[i] + 1)


class Parameters:

    def __init__(self, args: List[str]):
        # default parameter values
        self.approach = "branch-and-adjust"
        self.branching = "probabilities"
        self.delta = 0.001
        self.emphasis = 3
        self.instance_file = None
        self.output = "output"
        self.solution = ""
        self.threads = 1
        self.timelimit = None

        if not self.read_args(args):
            self.print_usage()
            sys.exit(0)

    def print_params(self):
        print("Parameter values:")
        print("     approach.:", self.approach)
        print("     branching:", self.branching)
        print("     delta....:", self.delta)
        print("     emphasis.:", self.emphasis)
        print("     output...:", self.output)
        print("     solution.:", self.solution)
        print("     threads..:", self.threads)
        print("     timelimit:", self.timelimit)
        print()

    def print_usage(self):
        print("Usage: python3 wta_cplex.py <input> [options]")
        print("    <input> : Path of the problem input file.")
        print("")
        print("Options:")
        print("    -approach <approach>   : selected approach/formulation to execute; possible values are:")
        print("                             {branch-and-adjust, probchain, underapprox, upperapprox}")
        print("                             (default: branch-and-adjust).")
        print("    -branching <branching> : branching strategy from {cplex, probabilities} (default: probabilities).")
        print("    -delta <delta>         : delta value (default: 0.001).")
        print("    -emphasis <emphasis>   : cplex search emphasis (default: 3).")
        print("    -output <output>       : output folder to save the log and additional info (default: output).")
        print("    -solution <file>       : file to save final solution (default: null).")
        print("    -threads <threads>     : maximum number of threads (default: 1).")
        print("    -timelimit <timelimit> : runtime limit (default: inf).")
        print("")

    def read_args(self, args) -> bool:
        print("Arguments: %s\n" % " ".join(args))

        i = 1
        while i < len(sys.argv):
            if sys.argv[i] == "-approach":
                self.approach = sys.argv[i + 1]
                print("Approach value set to '%s'" % self.approach)
                i += 2
            elif sys.argv[i] == "-branching":
                self.branching = sys.argv[i + 1]
                print("Banching strategy set to %s" % self.branching)
                i += 2
            elif sys.argv[i] == "-delta":
                self.delta = float(sys.argv[i + 1])
                print("Delta value set to %f" % self.delta)
                i += 2
            elif sys.argv[i] == "-emphasis":
                self.emphasis = int(sys.argv[i + 1])
                print("Cplex emphasis set to %d" % self.emphasis)
                i += 2
            elif sys.argv[i] == "-output":
                self.output = sys.argv[i + 1]
                print("Output folder set to %s" % self.output)
                i += 2
            elif sys.argv[i] == "-solution":
                self.solution = sys.argv[i + 1]
                print("Output solution file set to %s" % self.solution)
                i += 2
            elif sys.argv[i] == "-threads":
                self.threads = int(sys.argv[i + 1])
                print("Number of threads set to %d" % self.threads)
                i += 2
            elif sys.argv[i] == "-timelimit":
                self.timelimit = int(sys.argv[i + 1])
                print("Runtime limit set to %d seconds" % self.timelimit)
                i += 2
            elif not self.instance_file:
                self.instance_file = sys.argv[i]
                print("Reading input file %s" % self.instance_file)
                i += 1
            else:
                print("\nERROR: Unrecognized argument '%s'\n" % sys.argv[i], file=sys.stderr)
                return False

        if not self.instance_file:
            print("ERROR: No instance file provided!\n", file=sys.stderr)
            return False

        print()
        return True


class WTA:

    def compute_b(self, delta: float):
        self.b = {}
        self.B = []

        if delta < EPS:
            return

        sum_t = 0
        for j in self.inst.T:
            # adding smallest breakpoint for target j
            t = 0
            prod = 1
            for i in self.inst.W:
                prod *= (1 - self.inst.p[i, j]) ** self.inst.mu[i]
            if prod <= 0:
                prod = MIN_EXP
                print("    using prod = %f for target %d" % (prod, j))
            b_t = log(prod)
            self.b[j, t] = b_t
            t += 1

            # adding other breakpoints till b_t==0
            while b_t < 0:
                step_1 = lambda x: exp(x) - (exp(b_t) + exp(b_t) * (x - b_t)) - delta
                # x = fsolve(step_1, 1)[0]
                x = float(broyden1(step_1, 1, f_tol=1e-15))
                b_t = x if x < 0 else 0
                self.b[j, t] = b_t
                t += 1

                if b_t < 0:
                    step_2 = lambda x: exp(x) - (exp(b_t) + exp(x) * (x - b_t)) + delta
                    # x = fsolve(step_2, 1)[0]
                    x = float(broyden1(step_2, 1, f_tol=1e-15))
                    b_t = x if x < 0 else 0
                    if b_t == 0:
                        self.b[j, t] = b_t
                        t += 1

                # copying breakpoint from previous target (if they are iddentical)
                if t == 3:
                    for j_ in range(0, j):
                        if self.b[j_, t - 2] == self.b[j, t - 2] and self.b[j_, t - 1] == self.b[j, t - 1]:
                            # print("shortcut for j = %d" % j)
                            for _ in self.B[j_][t:]:
                                b_t = self.b[j_, t]
                                self.b[j, t] = b_t
                                t += 1
                            assert (b_t == 0)
                            break

            self.B.append(list(range(t)))
            sum_t += t

        print("    Average number of breakpoins: %.1f" % (sum_t / float(len(self.inst.T))))

    def get_lbda_obj(self, j: int, t: int):
        if t == 0:
            return self.inst.w[j] * exp(self.b[j, 0])
        elif t < self.B[j][-1]:
            return self.inst.w[j] * (exp(self.b[j, t]) - self.delta)
        else:
            return self.inst.w[j]


class WTA_ProbChain(WTA):
    """
    Formulation from Section 4.2, using Probability Chains
    """

    def __init__(self, inst: Instance, params: Parameters):
        self.delta = params.delta
        self.inst = inst
        self.model = cplex.Cplex()
        self.model.objective.set_sense(self.model.objective.sense.minimize)

        # setting cplex parameters
        self.model.parameters.read.datacheck.set(0)
        # self.model.parameters.emphasis.mip.set(params.emphasis)
        self.model.parameters.threads.set(params.threads)
        if params.timelimit:
            self.model.parameters.timelimit.set(params.timelimit)

        # computing b and B sets
        print("Computing b and B sets")
        self.b = {}
        self.B = []
        self.compute_b(self.delta)

        # adding variables
        print("Creating x, y and z variables")
        self.x = {
            (i, j, r): self.model.variables.add(types=[self.model.variables.type.binary], names=["x(%d,%d,%d)" % (i, j, r)])[0]
            for i in inst.W
            for j in inst.T
            for r in inst.R(i)
        }
        self.y = {
            (k, j, r): self.model.variables.add(names=["y(%d,%d,%d)" % (k, j, r)], lb=[0])[0]
            for k in inst.W[1:]
            for j in inst.T
            for r in inst.R(k)
        }
        self.z = {
            (i, j): self.model.variables.add(names=["z(%d,%d)" % (i, j)], lb=[0],
                                             obj=[inst.w[j] if i == inst.W[-1] else 0])[0]
            for i in inst.W
            for j in inst.T
        }

        # lbda variables are created only when delta > 0
        self.lbda = {}
        if self.delta > EPS:
            print("Creating lbda variables")
            self.lbda = {
                (j, t): self.model.variables.add(names=["lbda(%d,%d)" % (j, t)], lb=[0])[0]
                for j in inst.T
                for t in self.B[j]
            }

        # aliases to make code easier to read
        model = self.model
        x, z, y, lbda = self.x, self.z, self.y, self.lbda
        b, B = self.b, self.B
        delta = self.delta

        # adding constraints (4.17) -> c1
        print("Creating constraints (4.17)")
        for i in inst.W:
            expr = cplex.SparsePair(
                ind=[x[i, j, r] for j in inst.T for r in inst.R(i)],
                val=[r for j in inst.T for r in inst.R(i)]
            )
            model.linear_constraints.add(lin_expr=[expr], senses=["L"], rhs=[inst.mu[i]], names=["c1_%d" % i])

        # adding constraints (4.18) -> c2
        print("Creating constraints (4.18)")
        for i in inst.W:
            for j in inst.T:
                expr = cplex.SparsePair(
                    ind=[x[i, j, r] for r in inst.R(i)],
                    val=[1 for r in inst.R(i)]
                )
                model.linear_constraints.add(lin_expr=[expr], senses=["L"], rhs=[1], names=["c2_%d_%d" % (i, j)])

        # adding constraints (4.19) -> c3
        print("Creating constraints (4.19)")
        for j in inst.T:
            expr = cplex.SparsePair(
                ind=[z[0, j]] + [x[0, j, r] for r in inst.R(0)],
                val=[1] + [(1 - (1 - inst.p[0, j]) ** r) for r in inst.R(0)]
            )
            model.linear_constraints.add(lin_expr=[expr], senses=["E"], rhs=[1], names=["c3_%d" % j])

        # adding constraints (4.20) -> c4
        print("Creating constraints (4.20)")
        for k in inst.W[1:]:
            for j in inst.T:
                expr = cplex.SparsePair(
                    ind=[y[k, j, r] for r in inst.R(k)] + [z[k, j], z[k - 1, j]],
                    val=[1 for r in inst.R(k)] + [1, -1]
                )
                model.linear_constraints.add(lin_expr=[expr], senses=["E"], rhs=[0], names=["c4_%d_%d" % (k, j)])

        # adding constraints (4.21) -> c5
        print("Creating constraints (4.21)")
        for k in inst.W[1:]:
            for j in inst.T:
                for r in inst.R(k):
                    expr = cplex.SparsePair(
                        ind=[y[k, j, r], x[k, j, r]],
                        val=[1, -(1 - (1 - inst.p[k, j]) ** r)]
                    )
                    model.linear_constraints.add(lin_expr=[expr], senses=["L"], rhs=[0], names=["c5_%d_%d_%d" % (k, j, r)])

        # adding constraints (4.22) -> c6
        print("Creating constraints (4.22)")
        for k in inst.W[1:]:
            for j in inst.T:
                for r in inst.R(k):
                    expr = cplex.SparsePair(
                        ind=[y[k, j, r], z[k - 1, j]],
                        val=[1, -(1 - (1 - inst.p[k, j]) ** r)]
                    )
                    model.linear_constraints.add(lin_expr=[expr], senses=["L"], rhs=[0], names=["c6_%d_%d_%d" % (k, j, r)])

        if lbda:
            # adding constraints (5.10) -> sc1
            print("Creating constraints (5.10)")
            for j in inst.T:
                expr = cplex.SparsePair(
                    ind=[z[inst.W[-1], j]] + [lbda[j, 0]] + [lbda[j, t] for t in B[j][1:-1]] + [lbda[j, B[j][-1]]],
                    val=[-1] + [exp(b[j, 0])] + [exp(b[j, t]) - delta for t in B[j][1:-1]] + [1]
                )
                model.linear_constraints.add(lin_expr=[expr], senses=["L"], rhs=[0], names=["sc1_%d" % j])

            # adding constraints (5.11) -> sc2
            print("Creating constraints (5.11)")
            for j in inst.T:
                expr = cplex.SparsePair(
                    ind=[lbda[j, t] for t in B[j]] + [x[i, j, r] for i in inst.W for r in inst.R(i)],
                    val=[b[j, t] for t in B[j]] + [(-1) * log(1 - inst.p[i, j]) * r for i in inst.W for r in inst.R(i)]
                )
                model.linear_constraints.add(lin_expr=[expr], senses=["E"], rhs=[0], names=["sc2_%d" % j])

            # adding constraints (5.12) -> sc3
            print("Creating constraints (5.12)")
            for j in inst.T:
                expr = cplex.SparsePair(
                    ind=[lbda[j, t] for t in B[j]],
                    val=[1 for t in B[j]]
                )
                model.linear_constraints.add(lin_expr=[expr], senses=["E"], rhs=[1], names=["sc3_%d" % j])


class WTA_UnderApprox(WTA):
    """
    Formulation from Section 5, using Probability Chains with Linear Approximation From Below
    """

    def __init__(self, inst: Instance, params: Parameters):
        self.solution = []
        self.ub = float("inf")

        self.delta = params.delta
        self.inst = inst
        self.model = cplex.Cplex()
        self.model.objective.set_sense(self.model.objective.sense.minimize)

        # setting cplex parameters
        self.model.parameters.emphasis.mip.set(params.emphasis)
        # self.model.parameters.mip.strategy.heuristiceffort.set(0)
        # self.model.parameters.mip.strategy.heuristicfreq.set(-1)
        self.model.parameters.threads.set(params.threads)
        if params.timelimit:
            self.model.parameters.timelimit.set(params.timelimit)

        # computing b and B sets
        self.b = {}
        self.B = []
        self.compute_b(delta)

        # adding variables
        self.x = {
            (i, j, r): self.model.variables.add(types=[self.model.variables.type.binary], names=["x(%d,%d,%d)" % (i, j, r)])[0]
            for i in inst.W
            for j in inst.T
            for r in inst.R(i)
        }
        # self.y = {
        #     (k,j,r): self.model.variables.add(names=["y(%d,%d,%d)" % (k,j,r)], lb=[0])[0]
        #     for k in inst.W[1:]
        #     for j in inst.T
        #     for r in inst.R(k)
        # }
        # self.z = {
        #     (i,j): self.model.variables.add(names=["z(%d,%d)" % (i,j)], obj=[inst.w[j]], lb=[0])[0]
        #     for i in inst.W
        #     for j in inst.T
        # }
        self.lbda = {
            (j, t): self.model.variables.add(names=["lbda(%d,%d)" % (j, t)], lb=[0], obj=[self.get_lbda_obj(j, t)])[0]
            for j in inst.T
            for t in self.B[j]
        }

        # aliases to make code easier to read
        model = self.model
        x, z, lbda = self.x, self.z, self.lbda
        b, B = self.b, self.B

        # adding constraints (4.17)
        for i in inst.W:
            expr = cplex.SparsePair(
                ind=[x[i, j, r] for j in inst.T for r in inst.R(i)],
                val=[r for j in inst.T for r in inst.R(i)]
            )
            model.linear_constraints.add(lin_expr=[expr], senses=["L"], rhs=[inst.mu[i]], names=["c1_%d" % i])

        # adding constraints (4.18)
        for i in inst.W:
            for j in inst.T:
                expr = cplex.SparsePair(
                    ind=[x[i, j, r] for r in inst.R(i)],
                    val=[1 for r in inst.R(i)]
                )
                model.linear_constraints.add(lin_expr=[expr], senses=["L"], rhs=[1], names=["c2_%d_%d" % (i, j)])

        # adding constraints (5.11)
        for j in inst.T:
            expr = cplex.SparsePair(
                ind=[lbda[j, t] for t in B[j]] + [x[i, j, r] for i in inst.W for r in inst.R(i)],
                val=[b[j, t] for t in B[j]] + [-log(1 - inst.p[i, j]) * r for i in inst.W for r in inst.R(i)]
            )
            model.linear_constraints.add(lin_expr=[expr], senses=["E"], rhs=[0], names=["sc2_%d" % j])

        # adding constraints (5.12)
        for j in inst.T:
            expr = cplex.SparsePair(
                ind=[lbda[j, t] for t in B[j]],
                val=[1 for t in B[j]]
            )
            model.linear_constraints.add(lin_expr=[expr], senses=["E"], rhs=[1], names=["sc3_%d" % j])


class WTA_BranchAdjust_Int(WTA):
    """
    Formulation from Section 3 with Linear Approximation From Below,
    i.e. using objective function from Section 5 - see Eq. (5.10)
    """

    class BranchCallback(cplex.callbacks.BranchCallback):
        def __call__(self):
            if self.get_node_data() == "incumbent!!!":
                new_branches = False
                self.set_node_data(None)

                objval = self.get_objective_value()
                # feas = self.get_feasibilities(0, len(self.wta.x)+1)
                # obj = self.get_objective_coefficients(0, len(self.wta.x)+1)
                lbs = self.get_lower_bounds(0, len(self.wta.x) + 1)
                ubs = self.get_upper_bounds(0, len(self.wta.x) + 1)
                x = self.get_values(0, len(self.wta.x) + 1)

                for v in self.sorted_vars:
                    # making sure the variable has integer value
                    assert (np.fabs(x[v] - np.floor(x[v])) < EPS or np.fabs(x[v] - np.ceil(x[v])) < EPS)

                    # if LB == UB, then there is nothing to do here...
                    if lbs[v] >= ubs[v] - EPS:
                        continue

                    val = np.round(x[v])

                    # if UB > x
                    if ubs[v] > val + EPS:
                        self.make_branch(objval, variables=[(v, "U", val)])  # var <= x
                        self.make_branch(objval, variables=[(v, "L", val + 1)])  # var >= x + 1
                        # print("Adding branches var[%d] <= %f and var[%d] >= %f" % (i, val, i, val+1))
                        new_branches = True
                        break

                    # if LB < x
                    elif lbs[v] < val - EPS:
                        self.make_branch(objval, variables=[(v, "U", val - 1)])  # var <= x - 1
                        self.make_branch(objval, variables=[(v, "L", val)])  # var >= x
                        # print("Adding branches var[%d] <= %f and var[%d] >= %f" % (i, val-1, i, val))
                        new_branches = True
                        break

                if not new_branches:
                    assert self.get_num_branches() == 0

    class BoundCallback(cplex.callbacks.HeuristicCallback):
        def __call__(self):
            if self.wta.solution:
                print("\nNew %s improving solution found with *correct* cost of %f (relaxed had cost %f)\n" %
                      (self.message, self.wta.ub, self.wta.ub_ua), flush=True)
                self.set_solution(self.wta.solution, objective_value=self.wta.ub)
                self.wta.solution = []

    class CorrectObjectiveCallback(cplex.callbacks.IncumbentCallback):
        def __call__(self):
            solution = self.get_values()
            obj = self.get_objective_value()
            obj_adjusted = self.correct_objective_value(solution)

            if obj_adjusted < obj:
                solution_source = "heuristic" if self.get_solution_source() == self.solution_source.heuristic_solution else "node"
                for out in [sys.stdout, sys.stderr]:
                    print("", file=out)
                    print("Error with instance %s (delta=%f)" % (self.wta.params.instance_file, self.wta.params.delta), file=out)
                    print("    Invalid objective value: %f < %f" % (obj_adjusted, obj), file=out)
                    print("    Solution obtained by '%s'" % solution_source, file=out)
                    print("", file=out)

                # try:
                #     sol_path = self.wta.params.instance_file + "." + str(self.count) +  ".solution"
                #     lp_path = self.wta.params.instance_file + "." + str(self.count) + ".lp"
                #     self.count += 1
                #
                #     self.wta.model.write(lp_path)
                #
                #     for out in [sys.stdout, sys.stderr]:
                #         print("", file=out)
                #         print("Error with instance %s" % self.wta.params.instance_file, file=out)
                #         print("    Invalid objective value: %f < %f" % (obj_adjusted, obj), file=out)
                #         print("    Solution written to %s" % sol_path, file=out)
                #         print("    Model written to %s" % lp_path, file=out)
                #         print("", file=out)
                #
                #     f = open(sol_path, "w")
                #     print("Instance = %s" % self.wta.params.instance_file, file=f)
                #     print("Delta = %f" % self.wta.params.delta, file=f)
                #     print("Objective value = %.15f" % obj, file=f)
                #     print("'True' objective value = %.15f" % obj_adjusted, file=f)
                #     print("\nSolution:", file=f)
                #     for i in self.wta.inst.W:
                #         for j in self.wta.inst.T:
                #             if solution[self.wta.x[i,j]]:
                #                 print("    x(%d,%d) = %.50f" % (i,j,solution[self.wta.x[i,j]]), file=f)
                #     print("\nLambdas (for checking purposes):", file=f)
                #     for j in self.wta.inst.T:
                #         for t in self.wta.B[j]:
                #             if solution[self.wta.lbda[j,t]]:
                #                 print("    lbda(%d,%d) = %.50f" % (j,t,solution[self.wta.lbda[j,t]]), file=f)
                #
                #     f.close()
                #
                # except:
                #     pass

            # updating best bound information
            if obj_adjusted < self.wta.ub:
                self.wta.ub = obj_adjusted
                self.wta.ub_ua = self.get_objective_value()
                self.wta.solution = (list(range(len(solution))), list(solution))
                self.wta.bound_callback.message = "heuristic" if self.get_solution_source() == self.solution_source.heuristic_solution \
                    else "node"

            # rejecting solution in case it was not produced by the HeuristicCallback
            if self.get_solution_source() != self.solution_source.user_solution:
                if self.get_solution_source() == self.solution_source.node_solution:
                    self.set_node_data("incumbent!!!")
                self.reject()

        def correct_objective_value(self, solution: List[float]):
            sum = 0
            for j in self.wta.inst.T:
                prod = 1
                for i in self.wta.inst.W:
                    prod *= (1 - self.wta.inst.p[i, j]) ** solution[self.wta.x[i, j]]
                sum += self.wta.inst.w[j] * prod
            return sum

    def __init__(self, inst: Instance, params: Parameters):
        self.solution = []
        self.ub = float("inf")
        self.ub_ua = float("inf")

        self.delta = params.delta
        self.inst = inst
        self.params = params
        self.model = cplex.Cplex()
        self.model.objective.set_sense(self.model.objective.sense.minimize)

        # setting cplex parameters
        self.model.parameters.emphasis.numerical.set(1)
        # self.model.parameters.read.datacheck.set(0)
        self.model.parameters.emphasis.mip.set(params.emphasis)
        # self.model.parameters.mip.strategy.heuristiceffort.set(0)
        # self.model.parameters.mip.strategy.heuristicfreq.set(-1)
        self.model.parameters.threads.set(params.threads)
        # self.model.parameters.mip.tolerances.mipgap.set(1e-6)
        if params.timelimit:
            self.model.parameters.timelimit.set(params.timelimit)

        # computing b and B sets
        print("Computing b and B sets")
        self.b = {}
        self.B = []
        self.compute_b(params.delta)

        # adding variables
        print("Creating variables")
        self.x = {
            (i, j): self.model.variables.add(names=["x(%d,%d)" % (i, j)], lb=[0], ub=[inst.mu[i]],
                                             types=[self.model.variables.type.integer])[0]
            for i in inst.W
            for j in inst.T
        }
        self.lbda = {
            (j, t): self.model.variables.add(names=["lbda(%d,%d)" % (j, t)], lb=[0],
                                             obj=[self.get_lbda_obj(j, t)])[0]
            for j in inst.T
            for t in self.B[j]
        }

        # aliases to make code easier to read
        model = self.model
        x, lbda = self.x, self.lbda
        b, B = self.b, self.B

        # adding constraints (3.21)
        print("Creating constraints (3.21)")
        for j in inst.T:
            expr = cplex.SparsePair(
                ind=[lbda[j, t] for t in B[j]] + [x[i, j] for i in inst.W],
                val=[b[j, t] for t in B[j]] + [log(1 - inst.p[i, j]) * (-1) for i in inst.W]
            )
            model.linear_constraints.add(lin_expr=[expr], senses=["E"], rhs=[0], names=["sc2_%d" % j])

        # adding constraints (3.22)
        print("Creating constraints (3.22)")
        for j in inst.T:
            expr = cplex.SparsePair(
                ind=[lbda[j, t] for t in B[j]],
                val=[1 for t in B[j]]
            )
            model.linear_constraints.add(lin_expr=[expr], senses=["E"], rhs=[1], names=["sc3_%d" % j])

        # adding constraints (3.23)
        print("Creating constraints (3.23)")
        for i in inst.W:
            expr = cplex.SparsePair(
                ind=[x[i, j] for j in inst.T],
                val=[1 for j in inst.T]
            )
            model.linear_constraints.add(lin_expr=[expr], senses=["L"], rhs=[inst.mu[i]], names=["c1_%d" % i])

        # setting callbacks
        self.obj_callback = self.model.register_callback(self.CorrectObjectiveCallback)
        self.obj_callback.wta = self
        self.obj_callback.count = 1
        self.bound_callback = self.model.register_callback(self.BoundCallback)
        self.bound_callback.wta = self
        if params.branching == "probabilities":
            self.branch_callback = self.model.register_callback(self.BranchCallback)
            self.branch_callback.wta = self
            self.branch_callback.sorted_vars = [v[0] for v in sorted([(self.x[i, j], i, j) for (i, j) in self.x],
                                                                     key=lambda k: -self.inst.p[k[1], k[2]])]

        print()

    def solve(self):
        start = self.model.get_time()
        self.model.solve()
        self.runtime = self.model.get_time() - start

    def write_solution(self, path: str):
        try:
            solution = self.model.solution.get_values()

            f = open(path, "w")
            print(f"Instance = {self.params.instance_file}", file=f)
            print(f"Delta = {self.params.delta}", file=f)
            print(f"Objective value = {self.model.solution.get_objective_value():.15f}", file=f)
            print(f"Total runtime = {self.runtime:.2f} secods", file=f)
            print("\nSolution:", file=f)
            for i in self.inst.W:
                for j in self.inst.T:
                    if solution[self.x[i, j]]:
                        print("    x(%d,%d) = %.1f" % (i, j, solution[self.x[i, j]]), file=f)
            print("\nLambdas (for checking purposes):", file=f)
            for j in self.inst.T:
                for t in self.B[j]:
                    if solution[self.lbda[j, t]]:
                        print("    lbda(%d,%d) = %.50f" % (j, t, solution[self.lbda[j, t]]), file=f)
            f.close()

        except:
            print(f"Error when saving solution file {path}")


def main():
    # default argument values
    params = Parameters(sys.argv)
    inst = Instance(params.instance_file)

    wta = None
    if params.approach == "branch-and-adjust":
        params.print_params()
        wta = WTA_BranchAdjust_Int(inst, params)
    elif params.approach == "probchain":
        wta = WTA_ProbChain(inst, params)
    elif params.approach == "underapprox":
        wta = WTA_UnderApprox(inst, params)
    else:
        print("ERROR: approach '%s' is not available", file=sys.stderr)

    # solving linear relaxation
    # print("")
    # relax = wta.model
    # relax.relax()
    # relax.optimize()
    # print("")
    # print("Linear relaxation value = %f" % relax.objective_bound)
    # print("")

    # lp = os.path.join(params.output, os.path.basename(params.instance_file) + ".lp.gz")
    # print("Writing LP file %s" % lp)
    # wta.model.write(lp)

    # solving integer problem
    print("", flush=True)
    wta.solve()
    print("")
    print("\nSolutions status:", wta.model.solution.get_status_string())
    print("Number of nodes :", wta.model.solution.progress.get_num_nodes_processed())
    print("\nFinal optimality gap    = %f" % wta.model.solution.MIP.get_mip_relative_gap())
    print("Final lower bound value = %f" % wta.model.solution.MIP.get_best_objective())
    print("Final upper bound value = %f" % wta.model.solution.get_objective_value())
    print("")

    if params.solution:
        wta.write_solution(params.solution)


if __name__ == '__main__':
    main()
