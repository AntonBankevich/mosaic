import datetime
import os
from StringIO import StringIO

from typing import BinaryIO, Iterable, Generator, Optional, Union, IO, TextIO

from common import basic


class TokenWriter:
    def __init__(self, handler, info = ""):
        # type: (Union[BinaryIO, StringIO], str) -> None
        self.handler = handler
        self.new_line = True
        self.info = info

    def writeToken(self, token):
        # type: (str) -> None
        assert token.find(" ") == -1 and token.find("\n") == -1 and token != "", "Token:" + token
        if not self.new_line:
            self.handler.write(" ")
        self.new_line = False
        self.handler.write(token)

    def writeInt(self, val):
        # type: (Union[int, float]) -> None
        if not self.new_line:
            self.handler.write(" ")
        self.handler.write(str(val))
        self.new_line = False

    def writeTokenLine(self, token):
        # type: (str) -> None
        if not self.new_line:
            self.handler.write(" ")
        self.new_line = False
        self.handler.write(str(token))
        self.newLine()

    def writeIntLine(self, val):
        # type: (int) -> None
        self.writeToken(str(val))

    def newLine(self):
        self.handler.write("\n")
        self.new_line = True

    def writeTokens(self, tokens):
        # type: (Iterable[str]) -> None
        self.writeToken("List")
        for token in tokens:
            self.writeToken(token)
        self.newLine()

class TokenReader:
    def __init__(self, handler):
        # type: (Union[BinaryIO, StringIO]) -> None
        self.handler = handler
        self.line = None
        self.pos = None

    def readToken(self):
        # type: () -> Optional[str]
        while self.line is None or self.pos == len(self.line):
            self.line = self.handler.readline().split()
            self.pos = 0
        self.pos += 1
        if self.line[self.pos - 1] == "None":
            return None
        else:
            return self.line[self.pos - 1]

    def readInt(self):
        # type: () -> int
        return int(self.readToken())

    def readTokens(self):
        # type: () -> Generator[str]
        check = self.readToken()
        assert check == "List"
        for token in self.line[self.pos:]:
            yield token
        self.pos = len(self.line)

class SaveHandler:
    def __init__(self, dir, clean = False):
        # type: (str, bool) -> None
        self.dir = dir
        if clean:
            basic.recreate(self.dir)
        else:
            basic.ensure_dir_existance(self.dir)
        self.cnt = 0
        for name in os.listdir(self.dir):
            num = basic.parseNumber(name)
            if num is not None and num >= self.cnt:
                self.cnt = num + 1

    def getWriter(self, suffix = None):
        name = str(self.cnt) + "_" + datetime.datetime.now().strftime('%Y.%m.%d.%H.%M')
        self.cnt += 1
        if suffix is not None:
            name += "_" + suffix
        fn = os.path.join(self.dir, name)
        return TokenWriter(open(fn, "w"), fn)