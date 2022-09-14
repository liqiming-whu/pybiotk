from __future__ import annotations

import itertools
import shlex
import subprocess
import sys
from collections import deque
from enum import Enum
from typing import (
    Set,
    List,
    Any,
    Generic,
    Iterable,
    Iterator,
    Tuple,
    TypeVar,
    Callable,
    Literal,
    Optional,
    overload
)

__all__ = [
    "stream_type",
    "stream",
]

_T = TypeVar("_T")
_S = TypeVar("_S")


class stream_type(Enum):
    CMD_IN = 0
    STD_IN = 1
    FILE_IN = 2


class _stream_of:
    @staticmethod
    def sys_in(cmd: str, shell: bool = False) -> Iterator[str]:
        """single cmd;no linux pip
        """
        cmd2 = None
        if shell:
            cmd2 = cmd
        else:
            cmd2 = shlex.split(cmd)
        with subprocess.Popen(cmd2, stdout=subprocess.PIPE, encoding='utf8', shell=shell) as proc:
            assert proc.stdout is not None
            for line in proc.stdout:
                line = line.rstrip("\n")
                if line:
                    yield line

    @staticmethod
    def std_in() -> Iterator[str]:
        for line in sys.stdin:
            line = line.rstrip("\n")
            if line:
                yield line

    @staticmethod
    def file_in(file: str) -> Iterator[str]:
        with open(file) as f:
            for line in f:
                line = line.rstrip("\n")
                if line:
                    yield line


def it_window(iterables: Iterable[_T], window: int = 2, step: int = 1) -> Iterator[Tuple[_T, ...]]:
    it = iter(iterables)
    d = deque(itertools.islice(it, window), maxlen=window)
    n = 0
    yield tuple(d)
    for x in it:
        d.append(x)
        n += 1
        if n == step:
            yield tuple(d)
            n = 0
    if n != 0 and step > n:
        for _ in range(step-n):
            d.popleft()
        if d:
            yield tuple(d)


def _groupby(it: Iterable[_T], fun: Optional[Callable[[_T, _T], bool]] = None, compare: Literal['start', 'end'] = 'start') -> Iterator[Tuple[_T, ...]]:
    # it must be sorted before
    assert fun is not None
    it = iter(it)
    x = deque(itertools.islice(it, 1))
    compared = 0 if compare == 'start' else -1
    for i in it:
        if fun(x[compared], i):
            x.append(i)
        else:
            yield tuple(x)
            x.clear()
            x.append(i)
    if x:
        yield tuple(x)


class stream(Generic[_T]):

    def __init__(self, iterables: Iterable[_T]):
        self.it = iterables
        self.is_consumed = False

    def filter(self, func: Callable[[_T], bool]):
        self.it = filter(func, self.it)
        return self

    def flatten(self):
        return self.transform(itertools.chain.from_iterable)

    def concat(self, it: Iterable[_T]):
        self.it = itertools.chain(self.it, it)
        return self

    def chain(self):
        return self.transform(itertools.chain)

    def head(self, n: int = 10):
        self.it = itertools.islice(self.it, n)
        return self

    def abandon(self, n: int):
        self.it = itertools.islice(self.it, n, None)
        return self

    def tail(self, n: int = 10):
        self.it = iter(deque(self.it, maxlen=n))
        return self

    def move(self) -> Iterable[_T]:
        return self.it

    def list(self) -> List[_T]:
        return list(self.move())

    def set(self) -> Set[_T]:
        return set(self.move())

    def tuple(self) -> Tuple[_T, ...]:
        return tuple(self.move())

    def collect(self, kind: Literal['list', 'dict', 'tuple', 'set']):
        return eval(kind)(self.move())

    def to_file(self, name: str,  **kargs):
        with open(name, 'wt', **kargs) as f:
            for i in self.move():
                f.write(f"{i}\n")

    def transform(self, func: Callable[[Iterable[_T]], Iterable[_S]], *args, **kargs) -> stream[_S]:
        return stream._transform(self, func, *args, **kargs)

    def map(self, func: Callable[[_T], _S]) -> stream[_S]:
        return stream._map(self, func)

    def round(self, func: Callable[..., bool], b: int = 1, a: int = 1):
        return self.window(b+a+1).filter(lambda x: func(x[b])).flatten()

    def peek(self, func: Callable[[_T], Any]):
        """debug
        """
        def inner(it: Iterable[_T]) -> Iterable[_T]:
            for i in it:
                try:
                    func(i)
                    yield i
                except Exception:
                    raise RuntimeError('{}'.format(i))
        self.it = inner(self.it)
        return self

    def print(self):
        self.for_each(print)

    def for_each(self, fun: Callable[[_T], Any]) :
        assert self.is_consumed is False
        self.is_consumed = True
        for i in self.it:
            fun(i)

    @overload
    def count(self, fun: Callable[[int], _S]) -> _S: ...

    @overload
    def count(self) -> int: ...

    def count(self, fun=None):
        value = sum(map(lambda x: 1, self.it))
        if fun is not None:
            return fun(value)
        else:
            return value

    def nothing(self, *args, **kargs):
        return self

    def window(self, window: int = 2, step: int = 1) -> stream[Tuple[_T, ...]]:
        assert window > 0 and step > 0
        return stream(it_window(self.it, window=window, step=step))

    def groupby(self, key=None):
        return self.transform(stream._groupby, key=key)

    def groupby2(self, fun: Callable[[_T, _T], bool], compare: Literal['start', 'end'] = 'start'):
        """
        Parameters:
        ----------
        compare: where the behind element compare to in one group

        Example:
        >>> a=stream([1,1,2,2,3,3,4,5,4,6,7,8]).groupby2(lambda x,y: y-x<2, compare='end').print()
        >>> (1, 1, 2, 2, 3, 3, 4, 5, 4)
        >>> (6, 7, 8)
        >>> a=stream([1,1,2,2,3,3,4,5,4,6,7,8]).groupby2(lambda x,y: y-x<2).print()
        >>> (1, 1, 2, 2)
        >>> (3, 3, 4)
        >>> (5, 4, 6)
        >>> (7, 8)
        """
        return self.transform(_groupby, fun=fun, compare=compare)

    def sort(self, key: Optional[Callable] = None, reverse: bool = False):
        """sort a small number of element

        Parameters
        ----------
        key : Callable, optional
            the value to be sorted, by default None
        reverse : bool, optional
            ascending , by default False
        """
        self.it: Iterable[_T] = iter(sorted(self.it, key=key, reverse=reverse))
        return self

    def uniq(self, key: Optional[Callable] = None):
        """apply to a sorted iter, remove repeat item.
        Args:
            key : Callable, Defaults to None.
                func that gets a object values.

        Examples:
        --------
        >>> a = [1, 1, 2, 3, 3 , 4, 4, 5, 4, 4]
        >>> stream(iter(a)).uniq().for_each(print)
        ... 1
        ... 2
        ... 3
        ... 4
        ... 5
        ... 4
        """
        self.it = (list(y)[0] for _, y in itertools.groupby(self.it, key=key))
        return self

    @staticmethod
    def _map(x: stream, func: Callable[[_T], _S]) -> stream[_S]:
        x.it = map(func, x.it)
        return x

    @staticmethod
    def _transform(x: stream, func: Callable[[Iterable[_T]], Iterable[_S]], *args, **kargs) -> stream[_S]:
        x.it = iter(func(x.it, *args, **kargs))
        return x

    @staticmethod
    def _groupby(x, key=None) -> Iterator[Tuple[_T, ...]]:
        it = (tuple(g) for _, g in itertools.groupby(x, key=key))
        return it

    @staticmethod
    def of(type: stream_type, s: Optional[str] = None, shell: bool = False) -> stream[str]:
        if type is stream_type.STD_IN:
            return stream(_stream_of.std_in())
        else:
            assert s is not None
            if type is stream_type.CMD_IN:
                return stream(_stream_of.sys_in(s, shell=shell))
            else:
                return stream(_stream_of.file_in(s))
