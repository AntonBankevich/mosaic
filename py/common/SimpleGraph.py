from typing import Dict, List, Generator

from common import basic, SeqIO
from common.disjoint_set import DisjointSet


def toint(s):
    if s.startswith("\""):
        s = s[1:]
    if s.endswith("\""):
        s = s[:-1]
    return int(s)

class Edge:
    def __init__(self, id, start, fin, len, label, seq):
        self.id = id # type: str
        self.seq = seq # type: str
        self.label = label # type: str
        self.start = start # type: str
        self.len = len # type: int
        self.end = fin # type: str
        self.unique = (self.label.find("black") != -1)
        self.cov = None


class Vertex:
    def __init__(self, id, label):
        self.id = id
        self.inc = [] # type: List[Edge]
        self.out = [] # type: List[Edge]
        self.label = label


class SimpleGraph:
    def __init__(self):
        self.v = dict() #type: Dict[str, Vertex]
        self.visited = set()
        self.e = dict() #type: Dict[str, Edge]

    def AddVertex(self, vid, label = ""):
        if vid not in self.v:
            self.v[vid] = Vertex(vid, label)
        return self.v[vid]

    def AddEdge(self, id, start, fin, len, label, seq = None):
        v1 = self.AddVertex(start)
        v2 = self.AddVertex(fin)
        edge = Edge(id, start, fin, len, label, seq)
        self.e[id] = edge
        v1.out.append(edge)
        v2.inc.append(edge)
        return edge


    def ReadDot(self, f):
        for s in open(f, "r").readlines():
            tmp = s
            if s[0] != "\"":
                continue
            s =s.strip().split()
            if len(s) < 2 or s[1] != "->":
                vid = toint(s[0])
                self.AddVertex(vid, tmp)
                continue
            v_from = toint(s[0])
            v_to = toint(s[2])
            label = tmp
            cov = basic.parseNumber(tmp, tmp.find("k "))
            l = int(float(tmp[tmp.find("\\l") + 2: tmp.find("k ")]) * 1000)
            id = tmp[tmp.find("id ") + 3: tmp.find("\\l")]
            edge = self.AddEdge(id, v_from, v_to, l, label)
            edge.cov = cov
        return self

    def Merge(self):
        todel = []
        for v in self.v.values():
            if len(v.inc) == 1 and len(v.out) == 1 and v.inc[0].id != v.out[0].id:
                edge1 = self.e[v.inc[0].id]
                edge2 = self.e[v.out[0].id]
                if edge1.seq is None or edge2.seq is None:
                    seq = None
                else:
                    seq = edge1.seq + edge2.seq
                newEdge = Edge(edge1.id + "," + edge2.id, edge1.start, edge2.end, edge1.len + edge2.len, edge1.label + edge2.label, seq)
                newEdge.cov = int(float(edge1.cov + edge1.len + edge2.cov * edge2.len) / newEdge.len)
                self.v[edge1.start].out.remove(edge1)
                self.v[edge2.end].inc.remove(edge2)
                del self.e[edge1.id]
                del self.e[edge2.id]
                self.e[newEdge.id] = newEdge
                self.v[newEdge.start].out.append(newEdge)
                self.v[newEdge.end].inc.append(newEdge)
                todel.append(v)
        for v in todel:
            del self.v[v.id]



    def FillSeq(self, f, numeric = True):
        for s in SeqIO.parse_fasta(open(f, "r")):
            if numeric:
                s.id = str(basic.parseNumber(s.id))
            if s.id in self.e:
                self.e[s.id].seq = s.seq
                self.e[s.id].len = len(s.seq)
            if "-" + s.id in self.e:
                self.e["-" + s.id].seq = basic.RC(s.seq)
                self.e["-" + s.id].len = len(s.seq)
        return self

    def ReadGFA(self, f):
        seqs = dict()
        v = DisjointSet()
        for s in open(f, "r").readlines():
            s = s.split()
            if s[0] == "S":
                seqs[s[1]] = s[2]
                v.add((True, s[1], True))
                v.add((True, s[1],False))
                v.add((False, s[1], True))
                v.add((False, s[1],False))
            elif s[0] == "L":
                v1 = (s[2] == "+", s[1], True)
                v2 = (s[4] == "+", s[3], False)
                v.union(v1, v2)
                v1, v2 = ((not v2[0], v2[1], not v2[2]), (not v1[0], v1[1], not v1[2]))
                v.union(v1, v2)
        ids = dict()
        cnt = 1
        for vid, vl in v.listComponenets():
            self.AddVertex(str(cnt), str(vid))
            ids[vid] = str(cnt)
            cnt += 1
        for eid, seq in seqs.items():
            self.AddEdge(eid, ids[v.get((True, eid, False))], ids[v.get((True, eid, True))], len(seq), eid, seq)
            self.AddEdge("-" + eid, ids[v.get((False, eid, False))], ids[v.get((False, eid, True))], len(seq), "-" + eid, basic.RC(seq))
        return self

    def Find(self, minlen, v):
        self.visited.add(v)
        for e in self.v[v].out:
            if e.len <= minlen and e.end not in self.visited:
                self.Find(minlen, e.end)
        for e in self.v[v].inc:
            if e.len <= minlen and e.start not in self.visited:
                self.Find(minlen, e.start)

    def Split(self, minlen):
        # type: (int) -> Generator[List[Vertex]]
        visited = set()
        for v in self.v.values():
            if v.id in visited:
                continue
            self.visited = set()
            self.Find(minlen, v.id)
            visited = visited.union(self.visited)
            yield list(self.visited)

    def Draw(self, comp, out):
        out.write("digraph {\n")
        out.write("nodesep = 0.5;\n")
        for vid in comp:
            v = self.v[vid]
            if v.label != "":
                out.write(v.label)
        for v in self.v.values():
            for e in v.out:
                if e.start in comp or e.end in comp:
                    out.write(e.label + "\n")
        out.write("}\n")
