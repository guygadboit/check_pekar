from pdb import set_trace as brk
import os


ROOT = "./cumulative_results/"
MIN_POLYTOMY = 100


class Clade:
	def __init__(self, size, subclade_sizes):
		self.size = size
		self.subclade_sizes = subclade_sizes

	def polytomy_size(self):
		return len(self.subclade_sizes)


class Cladogram:
	def __init__(self, run, num_muts, clades):
		self.run = run
		self.num_muts = num_muts
		self.clades = clades

	def num_leaves(self):
		return sum([c.size for c in self.clades])

	def polytomy_size(self):
		return len(self.clades)

	def biggest_clade(self):
		best, best_size = self.clades[0], 0
		for clade in self.clades:
			if clade.size > best_size:
				best, best_size = clade, clade.size
		return best


class Sim:
	"""Represents one run of their simulation"""
	def __init__(self, run, one_mut, two_muts):
		"""one_mut and two_muts are the 1-mut and 2-mut cladograms for a
		particular run of the simulation"""
		self.run = run
		self.one_mut = one_mut
		self.two_muts = two_muts


def read(fname):
	with open(fname) as fp:
		for line in fp:
			yield line.strip()


def parse(fname, num_muts):
	"""Return a Cladogram"""
	clades = []
	run = int(os.path.basename(fname).split('_')[0])
	for line in read(fname):
		size, subs = line.split(maxsplit=1)
		subs = subs.strip('[]').split(',')
		size = int(size)
		subs = [int(x.strip()) for x in subs]
		clades.append(Clade(size, subs))
	return Cladogram(run, num_muts, clades)


def parse_dir(d, num_muts):
	"""Return a list of Cladograms from a directory"""
	ret = []
	for fname in os.listdir(d):
		ret.append(parse("{}/{}".format(d, fname), num_muts))
	return ret


def make_sims(one_mut, two_muts):
	"""Given two list of Cladograms create a single list of Sims that pairs
	them up correctly."""
	ret = []

	one = {cg.run: cg for cg in one_mut}
	two = {cg.run: cg for cg in two_muts}

	# Might as well have our list of Sims sorted for sanity
	for k in sorted(one.keys()):
		ret.append(Sim(k, one[k], two[k]))

	return ret


def find_cc_counts(sims):
	"""Return the various things they count up to see if we get the same
	numbers (spoiler: we do)"""
	cc_count = 0
	cc_count_1perc = 0
	cc_count_30perc = 0
	cc_count_30perc_twoPolytomies = 0

	for sim in sims:
		cg = sim.one_mut

		if len(cg.clades) != 2:
			continue

		sizes = [c.size for c in cg.clades]
		total = sum(sizes)

		if any([s < 2 for s in sizes]):
			continue

		cc_count += 1

		if min(sizes) > 0.01 * total:
			cc_count_1perc += 1

		if min(sizes) > 0.3 * total:
			cc_count_30perc += 1
			if all([c.polytomy_size() >= MIN_POLYTOMY for c in cg.clades]):
				cc_count_30perc_twoPolytomies += 1

	return (cc_count, cc_count_1perc, cc_count_30perc,
			cc_count_30perc_twoPolytomies)


def find_ab_counts(sims):
	"""Return the various things they count up to see if we get the same
	numbers (we don't)"""
	ab_count_30perc = 0
	ab_count_30perc_polytomy = 0
	ab_count_30perc_twoPolytomies = 0

	for sim in sims:
		cg = sim.two_muts

		if len(cg.clades) == 0:
			continue

		# The 1-mut cladograms include everything. The 2-mut ones only include
		# the 2-mut clades. This is why the total leaf count and base polytomy
		# sizes come from the one-mut cladogram.
		num_leaves = sim.one_mut.num_leaves()
		base_polytomy_size = sim.one_mut.polytomy_size()

		biggest = cg.biggest_clade()
		if 0.3 * num_leaves <= biggest.size <= 0.7 * num_leaves:

			# Is the biggest clade between 30% and 70%?
			ab_count_30perc += 1

			# Does it also have a base polytomy?
			if base_polytomy_size >= MIN_POLYTOMY:
				ab_count_30perc_polytomy += 1

				# And does that biggest clade itself have a polytomy?
				if biggest.polytomy_size() >= MIN_POLYTOMY:
					ab_count_30perc_twoPolytomies += 1

	return (ab_count_30perc, ab_count_30perc_polytomy,
			ab_count_30perc_twoPolytomies)


def main():
	one_mut = parse_dir(ROOT + "clade_analyses_CC", 1)
	two_muts = parse_dir(ROOT + "clade_analyses_AB", 2)
	sims = make_sims(one_mut, two_muts)

	# They get 116, 74, 16, 0. So do we.
	cc_counts = find_cc_counts(sims)
	print("cc_counts", cc_counts)

	# They get 119, 45, 5. We get 119, 45, 76
	ab_counts = find_ab_counts(sims)
	print("ab_counts", ab_counts)


if __name__ == "__main__":
	main()
