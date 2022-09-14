#!/usr/bin/env python3
"""Module enabling a sh like infix syntax (using pipes).
"""
import builtins
import functools
import gzip
import itertools
import os
import shlex
import subprocess
import sys
from collections import deque
from io import TextIOWrapper
from typing import (
    Set,
    List,
    Deque,
    Dict,
    Any,
    Tuple,
    Generic,
    Iterable,
    Iterator,
    TypeVar,
    Callable,
    Literal,
    Optional,
    Union,
    TextIO,
    overload
)

__all__ = [
    "cat",
    "zcat",
    "stdin",
    "cmdin",
    "mkdir",
    "cd",
    "Pipe",
    "head",
    "tail",
    "skip",
    "filter",
    "dedup",
    "uniq",
    "count",
    "to_list",
    "to_tuple",
    "to_dict",
    "to_set",
    "permutations",
    "join",
    "stdout",
    "write",
    "tee",
    "mapwith",
    "transform",
    "apply",
    "for_each",
    "flatten",
    "take_while",
    "drop_while",
    "skip_while",
    "kgroupby",
    "groupby",
    "groupby2",
    "sort",
    "reverse",
    "transpose",
    "window",
    "round",
    "chain",
    "concat",
    "islice"
]

_T = TypeVar("_T")
_S = TypeVar("_S")


class Pipe(Generic[_T]):
    """
    Represent a Pipeable Element :
    Described as :
    first = Pipe(lambda iterable: next(iter(iterable)))
    and used as :
    print [1, 2, 3] | first
    printing 1
    Or represent a Pipeable Function :
    It's a function returning a Pipe
    Described as :
    select = Pipe(lambda iterable, pred: (pred(x) for x in iterable))
    and used as :
    print [1, 2, 3] | select(lambda x: x * 2)
    # 2, 4, 6
    """

    def __init__(self, func: Callable[[Iterable[_T]], Any]):
        self.func = func
        functools.update_wrapper(self, func)

    def __ror__(self, other: Iterable[_T]):
        return self.func(other)

    def __call__(self, *args, **kwargs):
        return Pipe(lambda x: self.func(x, *args, **kwargs))


def cat(file: str) -> Iterator[str]:
    with open(file) as f:
        for line in f:
            line = line.rstrip("\n")
            if line:
                yield line


def zcat(file: str) -> Iterator[str]:
    with gzip.open(file, "rb") as f:
        for line in f:
            line = line.decode().rstrip("\n")
            if line:
                yield line


def stdin() -> Iterator[str]:
    for line in sys.stdin:
        line = line.rstrip("\n")
        if line:
            yield line


def cmdin(cmd: str, shell: bool = False) -> Iterator[str]:
    if shell:
        cmd2 = cmd
    else:
        cmd2 = shlex.split(cmd)
    with subprocess.Popen(cmd2, stdout=subprocess.PIPE, encoding='utf8', shell=shell) as proc:
        if proc.stdout is None:
            raise RuntimeError(proc.stderr)
        for line in proc.stdout:
            line = line.rstrip("\n")
            if line:
                yield line


def mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def cd(path):
    if os.path.exists(path):
        os.chdir(path)
    else:
        raise RuntimeError(f'{path} not exist')


@Pipe
def head(iterable: Iterable[_T], n: int = 10) -> Iterator[_T]:
    """Yield n of elements in the given iterable."""
    return itertools.islice(iterable, n)


@Pipe
def tail(iterable: Iterable[_T], n: int = 10) -> Deque:
    """Yield n of elements in the given iterable."""
    return deque(iterable, maxlen=n)


@Pipe
def skip(iterable: Iterable[_T], n: int) -> Iterator[_T]:
    """Skip n elements in the given iterable, then yield others."""
    return itertools.islice(iterable, n, None)


@Pipe
def filter(iterable: Iterable[_T], func: Callable[[_T], bool]) -> Iterable[_T]:
    return builtins.filter(func, iterable)


@Pipe
def dedup(iterable: Iterable[_T], key: Callable = lambda x: x) -> Iterator[_T]:
    """Only yield unique items. Use a set to keep track of duplicate data."""
    seen = set()
    for item in iterable:
        dupkey = key(item)
        if dupkey not in seen:
            seen.add(dupkey)
            yield item


@Pipe
def uniq(iterable: Iterable[_T], key: Callable = lambda x: x) -> Iterator[_T]:
    """Deduplicate consecutive duplicate values."""
    iterator = iter(iterable)
    try:
        prev = next(iterator)
    except StopIteration:
        return
    yield prev
    prevkey = key(prev)
    for item in iterator:
        itemkey = key(item)
        if itemkey != prevkey:
            yield item
        prevkey = itemkey


@overload
def count(iterable: Iterable[_T]) -> int: ...


@overload
def count(iterable: Iterable[_T], func: Callable[[int], _S]) -> _S: ...


@Pipe
def count(iterable, func=None):
    """
    Iteration will be complete
    """
    value = sum(builtins.map(lambda _: 1, iterable))
    if func is not None:
        return func(value)
    else:
        return value


@Pipe
def to_list(iterable: Iterable[_T]) -> List[_T]:
    return list(iterable)


@Pipe
def to_tuple(iterable: Iterable[_T]) -> Tuple[_T, ...]:
    return tuple(iterable)


@Pipe
def to_dict(iterable: Iterable[_T]) -> Dict[_T, _T]:
    return dict(iterable)


@Pipe
def to_set(iterable: Iterable[_T]) -> Set[_T]:
    return set(iterable)


@Pipe
def permutations(iterable: Iterable[_T], r: Optional[int] = None):
    # permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
    # permutations(range(3)) --> 012 021 102 120 201 210
    for x in itertools.permutations(iterable, r):
        yield x


@Pipe
def join(iterable: Iterable[_T], separator: str = ", "):
    return separator.join(builtins.map(str, iterable))


@Pipe
def stdout(iterable: Iterable[_T], func: Callable[[_T], Any] = str, end: str = "\n"):
    for item in iterable:
        sys.stdout.write(func(item) + end)


@Pipe
def write(iterable: Iterable[_T], path_or_buf: Union[str, TextIO], func: Callable[[_T], Any] = str, end: str = "\n"):
    if isinstance(path_or_buf, TextIOWrapper):
        ostream = path_or_buf
    else:
        ostream = open(path_or_buf, "w")
    for item in iterable:
        ostream.write(func(item) + end)
    ostream.close()


@Pipe
def tee(iterable: Iterable[_T], path_or_buf: Union[str, TextIO] = os.devnull, func: Callable[[_T], Any] = str, end: str = "\n"):
    if isinstance(path_or_buf, TextIOWrapper):
        if path_or_buf.name == "<stdout>" or path_or_buf.name == "/dev/stdout":
            ostream = open(os.devnull, "w")
        else:
            ostream = path_or_buf
    else:
        if path_or_buf == "/dev/stdout":
            ostream = open(os.devnull, "w")
        else:
            ostream = open(path_or_buf, "w")
    for item in iterable:
        ostream.write(func(item) + end)
        sys.stdout.write(func(item) + end)
    ostream.close()


@Pipe
def mapwith(iterable: Iterable[_T], func: Callable[[_T], _S]) -> Iterable[_S]:
    return builtins.map(func, iterable)


@Pipe
def transform(iterable: Iterable[_T], func: Callable[[Iterable[_T]], Iterable[_S]], *args, **kargs) -> Iterator[_S]:
    return iter(func(iterable, *args, **kargs))


@Pipe
def apply(iterable: Iterable[_T], func: Callable[[_T], Any]) -> Iterator[_T]:
    def inner(_iterable: Iterable[_T]):
        for i in _iterable:
            try:
                func(i)
                yield i
            except Exception:
                raise RuntimeError(str(i))

    return inner(iterable)


@Pipe
def for_each(iterable: Iterable[_T], func: Callable[[_T], Any]):
    """
    Iteration will be complete
    """
    for i in iterable:
        func(i)


@Pipe
def flatten(iterable: Iterable[_T], ignore_types=(str, bytes)) -> Iterator[_S]:
    for x in iterable:
        if isinstance(x, Iterable) and not isinstance(x, ignore_types):
            yield from x | flatten
        else:
            yield x


@Pipe
def take_while(iterable: Iterable[_T], func: Callable[[_T], bool]) -> Iterator[_T]:
    return itertools.takewhile(func, iterable)


@Pipe
def drop_while(iterable: Iterable[_T], func: Callable[[_T], bool]) -> Iterator[_T]:
    return itertools.dropwhile(func, iterable)


@Pipe
def skip_while(iterable: Iterable[_T], func: Callable[[_T], bool]) -> Iterator[_T]:
    for item in iterable:
        if func(item):
            continue
        yield item


@Pipe
def kgroupby(iterable: Iterable[_T], key: Optional[Callable[[_T], _S]] = None) -> Iterable[Tuple[_T, Tuple[_T, ...]]]:
    return ((i, tuple(g)) for i, g in itertools.groupby(iterable, key=key))


@Pipe
def groupby(iterable: Iterable[_T], key: Optional[Callable[[_T], _S]] = None) -> Iterable[Tuple[_T, ...]]:
    return (tuple(g) for _, g in itertools.groupby(iterable, key=key))


@Pipe
def groupby2(iterable: Iterable[_T], func: Callable[[_T, _T], bool], compare: Literal['start', 'end'] = 'start') -> Iterator[Tuple[_T, ...]]:
    """
    iterable must be sorted before
    Parameters:
    ----------
    compare: where the behind element compare to in one group

    Example:
    [1, 1, 2, 2, 3, 3, 4, 5, 4, 6, 7, 8] | groupby2(lambda x, y: y-x<2, compare='end') | stdout
    (1, 1, 2, 2, 3, 3, 4, 5, 4)
    (6, 7, 8)
    [1, 1, 2, 2, 3, 3, 4, 5, 4, 6, 7, 8] | groupby2(lambda x, y: y-x<2) | stdout
    (1, 1, 2, 2)
    (3, 3, 4)
    (5, 4, 6)
    (7, 8)
    """
    assert func is not None
    iterable = iter(iterable)
    x = deque(itertools.islice(iterable, 1))
    compared = 0 if compare == 'start' else -1
    for i in iterable:
        if func(x[compared], i):
            x.append(i)
        else:
            yield tuple(x)
            x.clear()
            x.append(i)
    if x:
        yield tuple(x)


@Pipe
def sort(iterable: Iterable[_T], key: Optional[Callable[[_T], _S]] = None, reverse: bool = False) -> Iterable[_T]:
    return sorted(iterable, key=key, reverse=reverse)


@Pipe
def reverse(iterable: Iterable[_T]) -> Iterable[_T]:
    return reversed(iterable)


@Pipe
def transpose(iterable: Iterable[_T]) -> Iterable[_T]:
    return builtins.zip(*iterable)


@Pipe
def window(iterables: Iterable[_T], window: int = 2, step: int = 1) -> Iterator[Tuple[_T, ...]]:
    assert window > 0 and step > 0
    iterable = iter(iterables)
    d = deque(itertools.islice(iterable, window), maxlen=window)
    n = 0
    yield tuple(d)
    for x in iterable:
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


@Pipe
def round(iterable: Iterable[_T], func: Callable[..., bool], b: int = 1, a: int = 1) -> Iterator[_S]:
    return iterable | window(b+a+1) | filter(lambda x: func(x[b])) | flatten


chain = Pipe(itertools.chain.from_iterable)
concat = Pipe(itertools.chain)
islice = Pipe(itertools.islice)
