#!/usr/bin/env python3
import sys
from abc import ABC, abstractmethod
from typing import Dict, List, Literal

import numpy as np
import pyBigWig
from numpy.lib.stride_tricks import sliding_window_view


class TrackFile(ABC):

    @property
    @abstractmethod
    def dchroms(self) -> Dict[str, int]: ...

    @abstractmethod
    def close(self): ...

    @abstractmethod
    def stats(self, chrom: str, start: int, end: int, stat: Literal['mean', 'max', 'min', 'sum', 'coverage', 'std'] = 'mean') -> np.ndarray:
        """
        assert turn None and np.nan to 0
        Args:
            chrom (str, optional): Defaults to None.
            start (int, optional): Defaults to None.
            end (int, optional):   Defaults to None.
            stat (Literal[, optional):  Defaults to 'mean'.

        Returns:
            np.ndarray: shape: (n,)
        """
        pass

    @abstractmethod
    def values(self, chrom: str, start: int, end: int) -> np.ndarray:
        """
        assert turn None and np.nan to 0
        Returns:
            np.ndarray: [, (n,end-start)]
        """
        pass

    @abstractmethod
    def __len__(self) -> int: ...

    def coverage(self, chrom: str, start: int, end: int) -> np.ndarray:
        """
        Returns:
            np.ndarray: shape: [-1, ]
        """
        a = self.values(chrom, start, end)
        return np.count_nonzero(a, axis=a.ndim-1) / a.shape[-1]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        if any([exc_type, exc_val, exc_tb]):
            return False


class Openbwn(TrackFile):
    def __init__(self, files: List[str]):
        self.files = files
        self.bws = [pyBigWig.open(x) for x in files]
        self.headers: List[Dict[str, int]] = [x.header() for x in self.bws]
        self.chroms: List[Dict[str, int]] = [x.chroms() for x in self.bws]

    @property
    def dchroms(self) -> Dict[str, int]:
        return self.chroms[0]

    def stats(self, chrom: str, start: int, end: int, stat: Literal['mean', 'max', 'min', 'sum', 'coverage', 'std'] = 'mean') -> np.ndarray:
        try:
            a = np.nan_to_num(np.array([bw.stats(chrom, start, end, exact=True, type=stat)[0] for bw in self.bws]))
            a[a == None] = 0
            return np.abs(a.astype(np.float32))
        except RuntimeError:
            sys.stderr.write(f'Invalid interval bounds! {chrom}:{start}-{end}')
            raise

    def values(self, chrom: str, start: int = None, end: int = None) -> np.ndarray:
        try:
            a = np.nan_to_num(np.array([bw.values(chrom, start, end) for bw in self.bws]))
            a[a == None] = 0
            return np.abs(a.astype(np.float32))
        except RuntimeError:
            sys.stderr.write(f'Invalid interval bounds! {chrom}:{start}-{end}')
            raise

    @staticmethod
    def coverage_bins_np_values(a: np.ndarray, nbins: int):
        assert a.shape[-1] % nbins == 0, 'length should be divisible by nbins '
        try:
            a2 = a.reshape((-1, int(a.shape[-1] / nbins)))
            return (np.count_nonzero(a2, axis=1) / a2.shape[-1]).reshape((a.shape[0], -1))
        except ValueError:
            sys.stderr.write(f"length of array is not divisible by nbins: {nbins}\n")
            raise

    def coverage_bins(self, chrom: str, start: int, end: int, nbins: int) -> np.ndarray:
        a = self.values(chrom, start, end)
        return self.coverage_bins_np_values(a, nbins)

    @staticmethod
    def scale_region_values_np_values(a: np.ndarray, length: int = 1000, nbins: int = 10) -> np.ndarray:
        """
        Returns:
            np.ndarray: shape: [-1, length/nbins]
        """
        n_base = length // nbins
        window_len = a.shape[-1] // n_base
        if window_len == 0:
            return np.zeros((a.shape[0], n_base))
        a = a[:, :window_len*n_base]
        return a.reshape(-1, window_len).mean(axis=1).reshape(a.shape[0], -1)

    def scale_region_values(self, chrom: str, start: int, end: int, length: int = 1000, nbins: int = 10):
        a = self.values(chrom, start, end)
        return self.scale_region_values_np_values(a, length, nbins)

    @staticmethod
    def coverage_sliding_window_np_values(a: np.ndarray, nbins: int) -> np.ndarray:
        a2 = sliding_window_view(a, (a.shape[0], nbins))
        a3 = (a2 > 0).sum(axis=3) / nbins
        return a3.squeeze(axis=0).transpose()

    def coverage_sliding_window(self, chrom: str, start: int, end: int, nbins: int):
        a = self.values(chrom, start, end)
        return self.coverage_sliding_window_np_values(a, nbins)

    def close(self):
        [x.close() for x in self.bws]

    def __len__(self) -> int:
        return len(self.bws)


class Openbw(TrackFile):
    def __init__(self, file: str):
        self.file = file
        self.bw = pyBigWig.open(file)
        self.header: Dict[str, int] = self.bw.header()
        self.chroms: Dict[str, int] = self.bw.chroms()

    @property
    def dchroms(self) -> Dict[str, int]:
        return self.chroms

    def stats(self, chrom: str, start: int, end: int, stat: Literal['mean', 'max', 'min', 'sum', 'coverage', 'std'] = 'mean') -> np.ndarray:
        try:
            a = np.nan_to_num(np.array(self.bw.stats(chrom, start, end, exact=True, type=stat)[0]))
            if a is None:
                a = np.float32(0)
            return np.abs(a.astype(np.float32))
        except RuntimeError:
            sys.stderr.write(f'Invalid interval bounds! {chrom}:{start}-{end}')
            raise

    def values(self, chrom: str, start: int = None, end: int = None) -> np.ndarray:
        try:
            a = np.nan_to_num(np.array(self.bw.values(chrom, start, end)))
            a[a == None] = 0
            return np.abs(a.astype(np.float32))
        except RuntimeError:
            sys.stderr.write(f'Invalid interval bounds! {chrom}:{start}-{end}')
            raise

    @staticmethod
    def coverage_bins_np_values(a: np.ndarray, nbins: int):
        assert a.shape[-1] % nbins == 0, 'length should be divisible by nbins '
        try:
            a2 = a.reshape((-1, int(a.shape[-1] / nbins)))
            return np.count_nonzero(a2, axis=1) / a2.shape[-1]
        except ValueError:
            sys.stderr.write(f"length of array is not divisible by nbins: {nbins}\n")
            raise

    def coverage_bins(self, chrom: str, start: int, end: int, nbins: int) -> np.ndarray:
        a = self.values(chrom, start, end)
        return self.coverage_bins_np_values(a, nbins)

    @staticmethod
    def scale_region_values_np_values(a: np.ndarray, length: int = 1000, nbins: int = 10) -> np.ndarray:
        """
        Returns:
            np.ndarray: shape: [-1, length/nbins]
        """
        n_base = length // nbins
        window_len = a.shape[-1] // n_base
        if window_len == 0:
            return np.zeros(n_base)
        a = a[:window_len*n_base]
        return a.reshape(-1, window_len).mean(axis=1)

    def scale_region_values(self, chrom: str, start: int, end: int, length: int = 1000, nbins: int = 10):
        a = self.values(chrom, start, end)
        return self.scale_region_values_np_values(a, length, nbins)

    @staticmethod
    def coverage_sliding_window_np_values(a: np.ndarray, nbins: int) -> np.ndarray:
        a2 = sliding_window_view(a, nbins)
        return np.count_nonzero(a2, axis=1) / nbins

    def coverage_sliding_window(self, chrom: str, start: int, end: int, nbins: int):
        a = self.values(chrom, start, end)
        return self.coverage_sliding_window_np_values(a, nbins)

    def close(self):
        self.bw.close()

    def __len__(self) -> int:
        return 1
