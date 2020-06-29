import re

from typing import Generator, Tuple, Optional

pattern = re.compile('([0-9]*)([MIDNSHP])')

def RCCigar(cigar):
    # type: (str) -> str
    res = []
    for n, c in list(pattern.findall(cigar))[::-1]:
        if n:
            res.append(n + c)
        else:
            res.append(c)
    return "".join(res)


def ReverseCigar(cigar):
    # type: (str) -> str
    res = []
    for c in cigar:
        if c == "D":
            c = "I"
        elif c == "I":
            c = "D"
        res.append(c)
    return "".join(res)


def CigarToList(cigar):
    # type: (str) -> Generator[Tuple[int, str]]
    if cigar == "=":
        yield (0, "=")
        return
    if cigar == "X":
        return
    for n, c in pattern.findall(cigar):
        if n:
            n = int(n)
        else:
            n = 1
        yield (n, c)

def CigarLen(cigar):
    res = 0
    for n, c in CigarToList(cigar):
        if c in "MD":
            res += n
    return res