import math

from typing import BinaryIO

from common import basic


class ScoringModel:
    def __init__(self):
        self. inf = 10000000000000

    def scoreIns(self, char):
        assert False

    def scoreDel(self, char):
        assert False

    def scoreMM(self, char1, char2):
        assert False

    def scoreHomo(self, char):
        assert False

    def scoreMatch(self, char):
        assert False

    def avgInsScore(self):
        assert False

    def avgDelScore(self):
        assert False


class SimpleScores(ScoringModel):
    def __init__(self, scoreIns = 0, scoreDel = 0, scoreMM = 0, scoreHomo = 0):
        ScoringModel.__init__(self)
        self.sIns = scoreIns
        self.sDel = scoreDel
        self.sMM = scoreMM
        self.sHomo = scoreHomo

    def scoreIns(self, char):
        return self.sIns

    def scoreDel(self, char):
        return self.sDel

    def scoreMM(self, char1, char2):
        return self.sMM

    def scoreHomo(self, char):
        return self.sHomo

    def scoreMatch(self, char):
        return 0

    def avgInsScore(self):
        return self.sIns

    def avgDelScore(self):
        return self.sDel


class ComplexScores(ScoringModel):
    def __init__(self):
        ScoringModel.__init__(self)
        self.sIns = [0.] * 4
        self.sDel = [0.] * 4
        self.sMM = [[0.]* 4 for i in range(4)]
        self.sHomo = [0.] * 4
        self.sMatch = [0.] * 4

    def load(self, handler):
        # type: (BinaryIO) -> None
        for s in handler.readlines():
            s = s.split()
            if len(s) == 0:
                continue
            if s[0] == "mat":
                self.sMatch[basic.letterToIndex(s[1])] = -math.log10(float(s[2]))
            elif s[0] == "mis":
                self.sMM[basic.letterToIndex(s[1][0])][basic.letterToIndex(s[1][-1])] = -math.log10(float(s[2]))
            elif s[0] == "del":
                self.sDel[basic.letterToIndex(s[1])] = -math.log10(float(s[2]))
                self.sHomo[basic.letterToIndex(s[1])] = -math.log10(float(s[2])) / 3
            elif s[0] == "ins":
                self.sIns[basic.letterToIndex(s[1])] = -math.log10(float(s[2]))

    def scoreIns(self, char):
        return self.sIns[basic.letterToIndex(char)]

    def avgInsScore(self):
        return sum(self.sIns) / 4

    def avgDelScore(self):
        return sum(self.sDel) / 4

    def scoreDel(self, char):
        return self.sDel[basic.letterToIndex(char)]

    def scoreMM(self, char1, char2):
        return self.sMM[basic.letterToIndex(char1)][basic.letterToIndex(char2)]

    def scoreHomo(self, char):
        return self.sHomo[basic.letterToIndex(char)]

    def scoreMatch(self, char):
        return self.sMatch[basic.letterToIndex(char)]