from typing import BinaryIO, Optional, Dict, Generator, Tuple, List

from common import basic

class EdgeInfo:
    def __init__(self, label, unique, cov, len, selfrc = False):
        # type: (str, bool, int, int, bool) -> EdgeInfo
        self.label = label
        self.unique = unique
        self.misc = []
        self.cov = cov
        self.selfrc = selfrc
        self.len = len

    def __repr__(self):
        return str([self.label, self.unique, self.cov, self.selfrc])

class DotParser:
    def __init__(self, dot):
        # type: (BinaryIO) -> DotParser
        self.dot = dot

    def parse(self, edge_ids = None):
        # type: (Optional[Dict[int, List[str]]]) -> Generator[Tuple[int, int, int, int, EdgeInfo]]
        for s in self.dot.readlines():
            if s.find("->") == -1:
                continue
            v_from = basic.parseNumber(s)
            v_to = basic.parseNumber(s, s.find("->"))
            eid = basic.parseNegativeNumber(s, s.find("id"))
            len = basic.parseNegativeNumber(s, s.find("\\l"))
            cov = basic.parseNumber(s, s.find("k "))
            unique = (s.find("black") != -1)
            src = (s.find("dir = both") != -1)
            if edge_ids is None or eid in edge_ids:
                if edge_ids is not None:
                    if "sink" in edge_ids[eid]:
                        v_to = "sink"
                    if "source" in edge_ids[eid]:
                        v_from = "source"
                yield eid, v_from, v_to, cov, EdgeInfo(s, unique, cov, len, src)