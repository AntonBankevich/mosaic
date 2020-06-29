from typing import Optional

from common import basic
from common.save_load import TokenWriter, TokenReader


class NamedSequence:
    def __init__(self, seq, id):
        # type: (str, str) -> NamedSequence
        self.id = id # type: str
        self.seq = seq.upper()

    def left(self):
        return 0

    def right(self):
        return len(self)

    def subSequence(self, left, right):
        # type: (int, int) -> NamedSequence
        return NamedSequence(self.seq[left:right], self.id + "(" + str(left) +"," + str(right) + ")")

    def __len__(self):
        # type: () -> int
        return len(self.seq)

    def RC(self):
        # type: () -> NamedSequence
        return NamedSequence(basic.RC(self.seq), basic.Reverse(self.id))

    def __getitem__(self, item):
        return self.seq[item]

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeToken(self.seq)
        handler.writeToken(self.id)

    @staticmethod
    def load(handler):
        # type: (TokenReader) -> NamedSequence
        return NamedSequence(handler.readToken(), handler.readToken())

    def __str__(self):
        return str(self.id) + "(" + str(self.__len__()) + ")"

    def __repr__(self):
        return str(self.id) + "(" + str(self.__len__()) + ")"

    def __eq__(self, other):
        return self.id == other.id

    def __neq__(self, other):
        return self.id != other.id


class SeqRecord(NamedSequence):
    def __init__(self, seq, id, qual = None, info = None):
        # type: (str, str, Optional[str], Optional[str]) -> SeqRecord
        assert qual is None or len(qual) == len(seq)
        NamedSequence.__init__(self, seq, id)
        self.qual = qual
        self.info = info

    def RC(self):
        # type: () -> SeqRecord
        qual = self.qual
        if qual is not None:
            qual = qual[::-1]
        return SeqRecord(basic.RC(self.seq), basic.Reverse(self.id), qual)

    def __getitem__(self, key):
        # type: (int) -> str
        return self.seq[key]

    def QualSubseq(self, l, r):
        # type: (int, int) -> Optional[str]
        if self.qual is not None:
            return self.qual[l: r]
        return None

    def subseq(self, l, r):
        # type: (int, int) -> SeqRecord
        if l != 0 or r != len(self.seq):
            return SeqRecord(self.seq[l:r], self.id + "(" + str(l + 1) +"-" + str(r) + ")", self.QualSubseq(l, r))
        else:
            return self