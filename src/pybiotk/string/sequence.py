# -*- coding: utf-8 -*-
import re

NOTATION = {
  "A": "A",
  "C": "C",
  "G": "G",
  "T": "T",
  "U": "U",
  "W": "[A|T|U]",
  "S": "[C|G]",
  "M": "[A|C]",
  "K": "[G|T|U]",
  "R": "[A|G]",
  "Y": "[C|T|U]",
  "B": "[C|G|T|U]",
  "D": "[A|G|T|U]",
  "H": "[A|C|T|U]",
  "V": "[A|C|G]",
  "N": "[A|C|G|T|U]",
  "Z": "0", 
}


def clean_string(string: str) -> str:
    return "".join(i for i in string.upper() if i in NOTATION)


def to_dna(string: str) -> str:
    return clean_string(string).replace("U", "T")


def to_rna(string: str) -> str:
    return clean_string(string).replace("T", "U")


def pattern(string: str) -> re.Pattern:
    return re.compile("".join(NOTATION.get(i, "") for i in string.upper()))
