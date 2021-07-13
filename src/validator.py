import sys


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
            w, t, p = int(splitted[0]), int(splitted[1]), float(splitted[2])
            self.p[w, t] = p

        # creating auxiliary sets
        self.W = list(range(self.m))
        self.T = list(range(self.n))


class Solution:

    def __init__(self, inst: Instance, path: str):
        self.instance_name = None
        self.delta = None
        self.objective = None

        self.inst = inst
        self.weapons = weapons = [[] for i in inst.W]
        self.targets = targets = [[] for j in inst.T]

        f = open(path, "r")

        i = 0
        lines = f.readlines()
        while i < len(lines):
            line = lines[i]

            # reading solution info
            if "Instance = " in line:
                self.instance_name = line.split()[-1]
            elif "Delta =" in line:
                self.delta = float(line.split()[-1])
            elif "Objective value = " in line:
                self.objective = float(line.split()[-1])

            # reading the solution itself
            elif "Solution:" in line:
                i += 1
                line = lines[i].strip()
                while i < len(lines) and line:
                    data = line.split()
                    w, t = data[0][2:-1].split(sep=",")
                    w, t = int(w), int(t)
                    value = int(round(float(data[-1])))
                    for v in range(value):
                        weapons[w].append(t)
                        targets[t].append(w)

                    i += 1
                    line = lines[i].strip()
                break
            i += 1

    def validate(self):
        # aliases (to make code easier to read)
        inst, weapons, targets = self.inst, self.weapons, self.targets

        # checking that weapons capacities are respected
        for w in inst.W:
            assert len(weapons[w]) <= inst.mu[w], \
                f"Error: capacity of weapon {w} is not respected ({len(weapons[w])} > {inst.mu[w]})"

        # computing survival probability of each target
        survival = []
        for t in inst.T:
            prod = 1
            for w in targets[t]:
                prod *= (1 - inst.p[w, t])
            survival.append(prod)
            print(f"Survival probability of target {t:3d}: {prod}")

        # computing overall solution cost
        cost = 0
        for t in inst.T:
            cost += inst.w[t] * survival[t]

        print("")
        print(f"Solution given cost...: {self.objective}")
        print(f"Solution computed cost: {cost}")
        print("")

        assert abs(cost - self.objective) < 1e-6, \
            f"Error: wrong solution cost, i.e. {self.objective} != {cost}"


def main():
    # default argument values
    instance_file = sys.argv[1]
    solution_file = sys.argv[2]
    inst = Instance(instance_file)
    sol = Solution(inst, solution_file)

    sol.validate()


if __name__ == '__main__':
    main()
