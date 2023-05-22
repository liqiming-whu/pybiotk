# -*- coding: utf-8 -*-
import re
from typing import (
    List,
    Literal,
    Optional,
    Iterator
)

from pybiotk.string.sequence import to_dna, to_rna, pattern


class Motif:
    def __init__(self, string: str, mtype: Literal['dna', 'rna'] = 'dna') -> None:
        self.string = string
        self.motif = to_dna(string) if mtype == 'dna' else to_rna(string)
        self.pattern = pattern(self.motif)
        
    def __repr__(self) -> str:
        return f"Motif({self.motif}):{self.pattern}"
    
    __str__ = __repr__
    
    def search(self, string: str) -> Optional[re.Match]:
        return self.pattern.search(string)
    
    def match(self, string: str) -> Optional[re.Match]:
        return self.pattern.match(string)
    
    def findall(self, string: str) -> List[str]:
        return self.pattern.findall(string)
    
    def finditer(self, string: str) -> Iterator:
        return self.pattern.finditer(string)
        