import itertools


class DisjointSet:
    def __init__(self):
        self.vals = dict()

    def add(self, item):
        self.vals[item] = item

    def get(self, item):
        if self.vals[item] != item:
            self.vals[item] = self.get(self.vals[item])
        return self.vals[item]

    def union(self, val1, val2):
        self.vals[self.get(val1)] = self.get(val2)

    def listComponenets(self):
        res = list(self.vals.keys())
        res = sorted(res, key = lambda item: self.get(item))
        return itertools.groupby(res, lambda item: self.get(item))

    def items(self):
        return list(self.vals.keys())