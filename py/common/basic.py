import heapq
import random
import shutil
import os
import sys
import time

from typing import Callable, Union, List

import common.log_params

rc = dict()
rc['A'] = 'T'
rc['C'] = 'G'
rc['G'] = 'C'
rc['T'] = 'A'
rc['a'] = 'T'
rc['c'] = 'G'
rc['g'] = 'C'
rc['t'] = 'A'
rc['N'] = 'N'


def RC(s):
    # type: (str) -> str
    res = list(s)[::-1]
    for i, c in enumerate(res):
        res[i] = rc[c]
    return "".join(res)


def ensure_dir_existance(dir):
    # type: (str) -> None
    if not os.path.exists(dir):
        os.makedirs(dir)


def recreate(dir):
    if os.path.exists(dir):
        shutil.rmtree(dir)
    ensure_dir_existance(dir)


def parseNumber(s, pos=0):
    # type: (str, int) -> Union[None, float, int]
    while pos < len(s) and s[pos] not in "0123456789":
        pos += 1
    if pos == len(s):
        return None
    pos1 = pos
    while pos1 < len(s) and s[pos1] in "0123456789.":
        pos1 += 1
    res = s[pos:pos1]
    if "." in res:
        return float(res)
    else:
        return int(res)


def parseNegativeNumber(s, pos=0):
    # type: (str, int) -> Union[None, float, int]
    while pos < len(s) and s[pos] not in "0123456789":
        pos += 1
    minus = False
    if pos > 0 and s[pos - 1] == '-':
        minus = True
    if pos == len(s):
        return None
    pos1 = pos
    while pos1 < len(s) and s[pos1] in "0123456789.":
        pos1 += 1
    if minus:
        pos -= 1
    res = s[pos:pos1]
    if "." in res:
        return float(res)
    else:
        return int(res)

def parseNegativeNumberAndMod(s, pos=0):
    # type: (str, int) -> str
    while pos < len(s) and s[pos] not in "0123456789":
        pos += 1
    minus = ""
    if pos > 0 and s[pos - 1] == '-':
        minus = "-"
    assert pos < len(s)
    return minus + s[pos:]

def best2(arr, better=lambda x, y: x < y):
    # type: (list[int], Callable[[int, int], bool]) -> tuple[int,int]
    assert len(arr) >= 2
    res = [0, 1]
    if better(arr[1], arr[0]):
        res = [1, 0]
    for i, val in enumerate(arr[2:], 2):
        if better(val, arr[res[1]]):
            res[1] = i
            if better(val, arr[res[0]]):
                res = res[::-1]
    return (res[0], res[1])


def merge(*iterables):
    prev = None
    for cur in heapq.merge(*iterables):
        if cur != prev:
            yield cur
            prev = cur


class OStreamWrapper:
    def __init__(self, *streams):
        self.streams = list(streams)
        self.active = True
        self.prefix = lambda x: ""
        self.level = common.log_params.LogPriority.log_level

    def fileno(self):
        return 1

    def write(self, string):
        if self.level >= common.log_params.LogPriority.common and self.active:
            for stream in self.streams:
                stream.write(string)

    def writelines(self, lines):
        if self.level >= common.log_params.LogPriority.common and self.active:
            for line in lines:
                for stream in self.streams:
                    stream.write(line)

    def info(self, *strings):
        self.log(common.log_params.LogPriority.main_parts, "INFO: " + " ".join(map(str, strings)))

    def trace(self, *strings):
        self.log(common.log_params.LogPriority.trace, "TRACE: " + " ".join(map(str, strings)))

    def substage(self, *strings):
        self.log(common.log_params.LogPriority.main_parts, "SUBSTAGE: " + " ".join(map(str, strings)))

    def warn(self, *strings):
        self.log(common.log_params.LogPriority.warning, "WARNING: " + " ".join(map(str, strings)))

    def log(self, level = 0, *strings):
        if level <= self.level and self.active:
            string = " ".join(map(str, strings))
            if self.active:
                self.write(self.prefix(string) + string + "\n")

    def block(self):
        self.active = False

    def release(self):
        self.active = True

    def flush(self):
        for stream in self.streams:
            stream.flush()


def Link(arr, dist):
    if len(arr) == 0:
        return
    arr = sorted([(pos, i) for i, pos in enumerate(arr)])
    left = arr[0]
    prev = arr[0]
    for val in arr[1:]:
        if val - prev > dist:
            yield left, prev
            left = val
        prev = val
    yield prev, arr[-1]


def Reverse(val):
    if isinstance(val, str):
        assert not val.startswith("--")
        if val.startswith("-"):
            return val[1:]
        else:
            return "-" + val
    elif isinstance(val, int):
        return -val
    assert False, "Tried to reverse an object that is neither number nor a string"


def Normalize(val):
    if isinstance(val, str):
        assert not val.startswith("--")
        if val.startswith("-"):
            return val[1:]
        else:
            return val
    elif isinstance(val, int):
        return abs(val)
    assert False, "Tried to normalize an object that is neither number nor a string"


def isCanonocal(val):
    if isinstance(val, str):
        assert not val.startswith("--")
        return not val.startswith("-")
    elif isinstance(val, int):
        return val > 0
    assert False, "Tried to check canonical an object that is neither number nor a string"

def canonical(val):
    if isCanonocal(val):
        return val
    else:
        return Reverse(val)


def quoted(val):
    return "\"" + str(val) + "\""


def fillRight(s, n, c=" "):
    while len(s) < n:
        s = s + c
    return s


def randomSequence(n):
    # type: (int) -> str
    res = ["A"] * n
    for i in range(n):
        res[i] = random.choice(["A", "C", "G", "T"])
    return "".join(res)


def parseLineName(name):
    # type: (str) -> List[str]
    if name[0] == "-":
        rc = True
        name = name[1:]
    else:
        rc = False
    if name.find("(") != -1:
        name = name[name.find("(") + 1: -1]
    seq = name.split(",")
    if rc:
        seq = map(Reverse, seq[::-1])
    return seq

def CreateLog(dir):
    old_logs_dir = os.path.join(dir, "old")
    ensure_dir_existance(old_logs_dir)
    log_file = os.path.join(dir, "log.info")
    num = len(filter(lambda name: name.endswith(".log"), os.listdir(old_logs_dir)))
    if os.path.isfile(log_file):
        shutil.copy(log_file, os.path.join(old_logs_dir, str(num) + ".log"))
    fasta_file = os.path.join(dir, "lines_knotted.fasta")
    if os.path.isfile(fasta_file):
        shutil.copy(fasta_file, os.path.join(old_logs_dir, str(num) + ".fasta"))
    log = open(log_file, "w")
    sys.stdout = OStreamWrapper(sys.stdout, log)
    sys.stdout.prefix = lambda s: time.strftime("%y.%m.%d %H:%M:%S") + "  "
    sys.stderr = sys.stdout

def letterToIndex(char):
    # type: (str) -> int
    if char == "A":
        return 0
    elif char == "C":
        return 1
    elif char == "G":
        return 2
    elif char == "T":
        return 3
    assert False
