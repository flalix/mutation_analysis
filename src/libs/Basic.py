#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
"""
# @Updated on 2022/12/06
# @Created on 2010/11/03
# @author: Flavio Lichtenstein
"""

import os  # sys

# import pickle5 as pickle
import pickle
import random
import re
import shutil
import zipfile  # , zlib
from collections import OrderedDict
from os.path import exists as exists
from os.path import join as osjoin
from pathlib import Path
from urllib.parse import urlparse
import requests

from typing import Any, List  # Optional, Iterable, Set, Tuple, ,

import numpy as np
import pandas as pd


class Basic(object):
    def __init__(self):
        self.log202 = np.log10(20) / np.log10(2)
        self.log2 = np.log10(2)
        self.ln2 = np.log(2)
        self.inv_ln2 = 1.0 / self.ln2
        self.ln2_2 = self.ln2**2
        self.convShannonTsallis = 1.417199516  # =1/LN(2)

        self.proteinMaxMER = np.log2(20)
        self.dnaMaxMER = np.log2(4)

    def getDnaNucleotides(self):
        return ["A", "T", "G", "C"]

    def getDnaNucleotideString(self):
        return "ATGC"

    def getRnaNucleotides(self):
        return ["A", "U", "G", "C"]

    def getrnaNucleotideString(self):
        return "AUGC"

    # in Hydropathy index[95] http://en.wikipedia.org/wiki/Amino_acid
    def getSeqAA(self):
        return [
            "A",
            "M",
            "C",
            "F",
            "L",
            "V",
            "I",
            "G",
            "T",
            "S",
            "W",
            "Y",
            "P",
            "H",
            "N",
            "D",
            "E",
            "Q",
            "K",
            "R",
        ]

    def getStringAA(self):
        return "AMCFLVIGTSWYPHNDEQKR"

    def getAaPos(self, string, aa):
        return string.find(aa)


class Log_writer(object):
    """log writer"""

    def __init__(self, fileLog="log.txt", path="../log/", remove=True):
        # global fileLog
        # fileLog = "../results/log.txt"
        self.fileLog = fileLog
        self.pathLog = path
        filename = osjoin(path, fileLog)

        if remove:
            if exists(filename):
                os.remove(filename)

    def write_log(self, text, to_append=True, verbose=False):  # encondig = 'utf-8'
        write_txt(text, self.fileLog, self.pathLog, to_append=to_append, verbose=verbose)


def simple_replace(stri):
    if stri is None:
        return None
    stri = stri.strip()
    if stri == "":
        return ""

    return stri.replace("  ", " ")


def download_url_file(url: str, fname:str, root_file: Path, 
                      timeout: int = 90, force: bool = False, verbose: bool = False) -> bool:
    root_file.mkdir(parents=True, exist_ok=True)

    filename = root_file / fname

    if filename.exists() and not force:
        if verbose:
            print(f"File already exists: {filename}")
        return True

    try:
        r = requests.get(url, stream=True, timeout=timeout)
        r.raise_for_status()
    except requests.RequestException as e:
        print(f"Error downloading image from {url}: {e}")
        return False

    try:
        with open(filename, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
    except Exception as e:
        print(f"Error saving image to {filename}: {e}")
        return False
    
    if verbose:
        print(f"Image downloaded successfully at {filename}.")

    return True

def prepare_id(_id):
    return _id.split("=")[0].strip().replace("/", "-").replace(" ", "_")


def prepare_title(title, maxCount=12):
    if title is None:
        return None
    title = title.strip()
    if title == "":
        return "???"

    title = title.replace("/", "-").replace("_", ": ")
    words = title.split(" ")

    stri = ""
    count = 0
    for word in words:
        stri += word + " "
        count += 1
        if count == maxCount:
            stri = stri.rstrip() + "\n"
            count = 0

    return stri.rstrip()


def full_replace(stri):
    if stri is None:
        return None
    if stri.strip() == "":
        return ""

    return (
        stri.replace("  ", "_")
        .replace(" ", "_")
        .replace("/", "-")
        .replace("[", "(")
        .replace("]", ")")
        .replace("\n", " - ")
    )


def full_replace_lower(stri):
    if stri is None:
        return None
    return full_replace(stri).lower()


def write_txt(text, filename, path=Path("./"), to_append=False, verbose=False) -> bool:
    if not path.exists():
        os.mkdir(path)

    filename = path / filename
    h = None
    try:
        ret = True
        if to_append:
            h = open(filename, mode="a+")  # encondig = 'utf-8',
        else:
            h = open(filename, mode="w")  # encondig = 'utf-8',
        h.write(text)
    except ValueError:
        print(f"Error '{ValueError}', writing in '{filename}'")
        ret = False
    finally:
        if h:
            h.close()

    if not ret:
        return False

    if verbose:
        print(f"File saved: {filename}")
    return True


"""
https://thispointer.com/5-different-ways-to-read-a-file-line-by-line-in-python/
"""


def read_txt(
    filename, path=Path("./"), iniLine=None, endLine=None, verbose=False
):  # encondig = 'utf-8',
    filename = path / filename
    text = []

    try:
        h = open(filename, mode="r")

        if iniLine is not None or endLine is not None:
            iniLine = -1 if iniLine is None else iniLine
            endLine = np.inf if endLine is None else endLine

            count = 0
            while True:
                if count < iniLine:
                    _ = h.readline()
                else:
                    line = h.readline()
                    if not line:
                        break

                    text.append(line.replace("\n", ""))

                count += 1
                if count > endLine:
                    break
        else:
            text = [line.replace("\n", "") for line in h.readlines()]

        h.close()
        if verbose:
            print("File read at '%s'" % (filename))
    except ValueError:
        print(f"Error '{ValueError}' while reading '{filename}'")

    return "\n".join(text)


# encondig = 'utf-8',
def pdreadcsv(
    fname: str,
    path: Path = Path("./"),
    sep: str = "\t",
    dtype: dict = {},
    colnames: List = [],
    index_col=-1,
    skiprows: int = 0,
    selcols: List = [],
    sortcols: List = [],
    low_memory: bool = False,
    removedup: bool = False,
    header: int = 0,
    verbose: bool = False,
) -> pd.DataFrame:

    if not path.exists():
        print(f"Path does not exists: '{path}'")
        return pd.DataFrame()

    filename = path / fname

    if not filename.exists():
        print(f"File does not exist: '{filename}'")
        return pd.DataFrame()

    try:
        if dtype == {}:
            if index_col < 0:
                df = pd.read_csv(
                    filename, sep=sep, low_memory=low_memory, skiprows=skiprows, header=header,
                )
            else:
                df = pd.read_csv(
                    filename, sep=sep, low_memory=low_memory, skiprows=skiprows, header=header, index_col=index_col,
                )
        else:
            if index_col < 0:
                df = pd.read_csv(
                    filename,
                    sep=sep,
                    dtype=dtype,
                    low_memory=low_memory,
                    skiprows=skiprows,
                    header=header,
                )
            else:
                df = pd.read_csv(
                    filename,
                    sep=sep,
                    dtype=dtype,
                    low_memory=low_memory,
                    skiprows=skiprows,
                    header=header,
                    index_col=index_col
                )
    except Exception as e:
        print(f"Error reading csv/tsv '{filename}': {e}")
        return pd.DataFrame()

    try:
        if len(colnames) > 0:
            df.columns = colnames
        if len(selcols) > 0:
            df = df[selcols]
        if len(sortcols) > 0:
            df = df.sort_values(sortcols)
        if removedup:
            df = df.drop_duplicates()
    except ValueError:
        print(f"Error '{ValueError}' columns/selecting/sorting '{filename}'")
        return df

    if verbose:
        print("Table opened (%s) at '%s'" % (df.shape, filename))
    return df


def month_to_num(stri):
    if isinstance(stri, int):
        return stri
    if not isinstance(stri, str):
        return stri

    stri2 = stri.lower()

    if stri2 == "jan":
        return 1
    if stri2 == "feb" or stri2 == "fev":
        return 2
    if stri2 == "mar":
        return 3
    if stri2 == "apr" or stri2 == "abr":
        return 4
    if stri2 == "may" or stri2 == "mai":
        return 5
    if stri2 == "jun":
        return 6
    if stri2 == "jul":
        return 7
    if stri2 == "aug" or stri2 == "ago":
        return 8
    if stri2 == "sep" or stri2 == "set":
        return 9
    if stri2 == "oct" or stri2 == "out":
        return 10
    if stri2 == "nov":
        return 11
    if stri2 == "dec" or stri2 == "dez":
        return 12

    return stri


def replace_space(stri, default: str = "_"):
    if not isinstance(stri, str):
        return stri
    return stri.replace(" ", default)


def remove_spaces(stri: str, replace_nans: bool = True):
    if replace_nans:
        stri = stri.replace("nan", "")
    stri = re.sub(" +", " ", stri)
    return stri.strip()


# encondig = 'utf-8',
def pdwritecsv(df, filename, path: Path = Path("./"), sep="\t", index=False, verbose=False):

    if not path.exists():
        print("Path does not exists: '%s'" % (path))
        return False

    filename = path / filename
    try:
        df.to_csv(filename, sep=sep, index=index)  # , encondig = encondig
    except ValueError:
        print(f"Error '{ValueError}', writing in '{filename}'")
        return False

    if verbose:
        print("Table saved (%s) at '%s'" % (df.shape, filename))
    return True


def dumpdic(dic, filename, path: Path = Path("./"), verbose=True):
    return pddumpdic(dic=dic, filename=filename, path=path, verbose=verbose)


def pddumpdic(dic, filename, path: Path = Path("./"), verbose=True):

    if not path.exists():
        print(f"Path does not exists: '{path}'")
        return False

    filename = path / filename
    try:
        # Store data (serialize)
        with open(filename, "wb") as handle:
            pickle.dump(dic, handle)  # , protocol=pickle.HIGHEST_PROTOCOL
    except ValueError:
        print(f"Error '{ValueError}' in pickle dump '{filename}'")
        return False

    if verbose:
        print("Dictionary saved at '%s'" % (filename))
    return True


def loaddic(filename, path="./", verbose=True):
    return pdloaddic(filename=filename, path=path, verbose=verbose)


def pdloaddic(filename: str, path: str = "./", verbose: bool = False) -> dict:

    dic = {}

    if not exists(path):
        print("Path does not exists: '%s'" % (path))
        return dic

    filename = osjoin(path, filename)

    if not exists(filename):
        print("File does not exist: '%s'" % (filename))
        return dic

    try:
        # Load data (deserialize)
        with open(filename, "rb") as handle:
            dic = pickle.load(handle)

        if verbose:
            print("Dictionary read at '%s'" % (filename))
    except ValueError as e:
        print(f"Error '{ValueError}' in pickle loading '{filename}': {e}")

    return dic


def padl(val, num, charac="0"):
    return pad(val, num, charac)


def pad(val, num, charac="0"):
    sval = str(val).strip()
    dif = num - len(sval)
    if dif <= 0:
        return sval

    zeros = [charac] * dif
    return "".join(zeros) + sval


def padr(val, num, charac="0"):
    sval = str(val).strip()
    dif = num - len(sval)
    if dif <= 0:
        return sval

    zeros = [charac] * dif
    return sval + "".join(zeros)


def try_float(num, except_ret=None):
    try:
        return float(num)
    except ValueError:
        return except_ret


def isfloat(num):
    try:
        num = float(num)
        return True
    except ValueError:
        return False


def try_int(num, except_ret=None):
    try:
        return int(num)
    except ValueError:
        return except_ret


def isint(num):
    try:
        num = int(num)
        return True
    except ValueError:
        return False


def isint_v2(x):
    try:
        a = float(x)
        b = int(a)
    except ValueError:
        return False

    return a == b


def return_integers(vals, except_ret=None):
    return [int(x) for x in vals if try_int(x, except_ret) is not None]


def return_floats(vals, except_ret=None):
    return [float(x) for x in vals if try_float(x, except_ret) is not None]


##--- only strings??
def merge_by_columns_inner_outer(dflista, keys, how, fillna=False):
    dfn = dflista[0]

    for i in range(1, len(dflista)):
        df2 = dflista[i]
        dfn = pd.merge(dfn, df2, how=how, on=keys)

    if fillna:
        dfn = dfn.fillna(0)

    return dfn


def merge_by_columns_inner(dflista, keys, fillna=False):
    dfn = dflista[0]

    for i in range(1, len(dflista)):
        df2 = dflista[i]
        dfn = pd.merge(dfn, df2, how="inner", on=keys)

    if fillna:
        dfn = dfn.fillna(0)

    return dfn


def merge_by_columns_outer(dflista, keys, fillna=False):
    return merge_by_columns(dflista, keys, fillna)


def merge_by_columns(dflista, keys, fillna=False):
    dfn = dflista[0]

    for i in range(1, len(dflista)):
        df2 = dflista[i]
        dfn = pd.merge(dfn, df2, how="outer", on=keys)

    if fillna:
        dfn = dfn.fillna(0)

    return dfn


def prepare_figname(figname, endtype="png"):
    posi = [i for i in range(len(figname)) if figname[i] == "."]
    if len(posi) > 0:
        figname = figname[: posi[0]]

    figname = title_replace(figname)

    figname = figname + "." + endtype
    return figname


def title_replace(title: str) -> str:
    title = title.strip()
    if title[-1] == ".":
        title = title[:-1]

    return (
        title.replace(",", "_")
        .replace(";", "_")
        .replace(" - ", "_")
        .replace("\n", "_")
        .replace("<br>", "_")
        .replace("  ", "_")
        .replace(" ", "_")
        .replace("/", "-")
        .replace("_=_", "_")
        .replace("=", "_")
        .replace("___", "_")
        .replace("__", "_")
        .replace(":", "")
    )


# is_in(['aaaa', 'bbbbb', 'abcd'], ['a','x'])
def is_in(series, lista):
    series = [x.lower().replace("-", "") for x in series]
    return [np.sum([1 if x in k else 0 for x in lista]) > 0 for k in series]


colorscale = [
    [1.0, "rgb(165,0,38)"],
    [0.8888888888888888, "rgb(215,48,39)"],
    [0.7777777777777778, "rgb(244,109,67)"],
    [0.6666666666666666, "rgb(253,174,97)"],
    [0.5555555555555556, "rgb(254,224,144)"],
    [0.4444444444444444, "rgb(224,243,248)"],
    [0.3333333333333333, "rgb(171,217,233)"],
    [0.2222222222222222, "rgb(116,173,209)"],
    [0.1111111111111111, "rgb(69,117,180)"],
    [0.0, "rgb(49,54,149)"],
]


def set_color_scale(min, max):
    maxi = max
    if abs(min) > maxi:
        maxi = abs(min)

    pos = np.linspace(-maxi, max, num=10)

    colorscale = [
        [pos[9], "rgb(165,0,38)"],
        [pos[8], "rgb(215,48,39)"],
        [pos[7], "rgb(244,109,67)"],
        [pos[6], "rgb(253,174,97)"],
        [pos[5], "rgb(254,224,144)"],
        [pos[4], "rgb(224,243,248)"],
        [pos[3], "rgb(171,217,233)"],
        [pos[2], "rgb(116,173,209)"],
        [pos[1], "rgb(69,117,180)"],
        [pos[0], "rgb(49,54,149)"],
    ]

    return colorscale


def columns_to_case(df, pre_cols, cols):

    dff = None
    for col in cols:
        dfa = df[pre_cols + [col]].copy()
        dfa["class"] = col
        dfa.columns = pre_cols + ["val", "class"]

        if dff is None:
            dff = dfa
        else:
            dff = dff.append(dfa)

    return dff


def char_frequency(mat: List):
    dic = OrderedDict()
    for c in mat:
        if c in dic.keys():
            dic[c] += 1
        else:
            dic[c] = 1
    return dic


def best_nucleotide(mat, with_gaps=False):
    dicc = char_frequency(mat)

    nuc_list = ["A", "T", "G", "C"]
    if with_gaps:
        nuc_list += ["-"]

    maxi = 0
    best_key = None
    for key, value in dicc.items():
        if key in nuc_list:
            if value > maxi:
                maxi = value
                best_key = key

    return best_key


def best_amino_acid(mat, with_gaps=False):
    dicc = char_frequency(mat)

    aa_list = [
        "A",
        "M",
        "C",
        "F",
        "L",
        "V",
        "I",
        "G",
        "T",
        "S",
        "W",
        "Y",
        "P",
        "H",
        "N",
        "D",
        "E",
        "Q",
        "K",
        "R",
    ]
    if with_gaps:
        aa_list += ["-"]

    maxi = 0
    best_key = None
    for key, value in dicc.items():
        if key in aa_list:
            if value > maxi:
                maxi = value
                best_key = key

    return best_key


import multiprocessing as mp


def get_cpus(decrement_cpus=1, verbose=True):
    maxcpus = mp.cpu_count()
    cpus = maxcpus - decrement_cpus

    if verbose:
        print("max CPUs %d - using %d" % (maxcpus, cpus))
    return cpus


def break_lines_length(
    stri: str, split_char: str = " ", sep: str = "<br>", maxLen: int = 30
) -> str:
    return break_line_per_length(stri, split_char, sep, maxLen)


def break_line_per_length(
    stri: str, split_char: str = " ", sep: str = "<br>", maxLen: int = 30
) -> str:
    mat = stri.split(split_char)
    text = ""
    temp = ""
    insert_break = False
    for term in mat:
        if insert_break:
            text += sep + term
            insert_break = False
        else:
            text += term

        temp += term

        if len(temp) > maxLen:
            temp = ""
            insert_break = True
        else:
            text += " "

    return text.strip()


def break_list(
    term_list: List, n_elems: int = 5, CR: str = "<br>", spacer: str = "; ", maxLines: int = -1
) -> str:

    iline = 0
    icount = 0
    text = ""
    for i in range(len(term_list)):
        term = term_list[i]
        icount += 1
        text += term
        if n_elems > 0 and icount == n_elems:
            icount = 0
            text += CR
            last_term_break = True
            iline += 1
        else:
            text += spacer
            last_term_break = False

        if maxLines > 0 and iline >= maxLines:
            print(">>> maxLines")
            text += "..."
            return text

    if last_term_break:
        text = text[: -len(CR)]
    else:
        text = text[: -len(spacer)]

    return text


def break_lines(
    stri: str,
    sep: str = " ",
    nwords: int = 25,
    CR: str = "<br>",
    maxLen: int = 250,
    maxLines: int = -1,
) -> str:
    stri = stri.strip()
    text = ""
    mat = stri.split(sep)
    lines = 0
    for i in range(len(mat), nwords):
        mat2 = mat[i : (i + nwords)]
        if len(mat2) == 0:
            break

        newtext = sep.join(mat2)
        if len(newtext) > maxLen:
            while True:
                text += newtext[:maxLen] + CR
                if len(newtext) <= maxLen:
                    break
                newtext = newtext[maxLen:]
                lines += 1
        else:
            text += newtext + CR
            lines += 1

        if maxLines > 0 and lines > maxLines:
            text += "..."
            return text

    return text[: -len(CR)]


def shuffle_nums(nrows, N, howManySets, verbose=False):

    samples = np.arange(0, nrows)
    random.shuffle(samples)

    # print(">> samples: %d x %d"%(nrows, len(samples)))

    sampList = []
    for i in range(howManySets):
        if nrows >= (i + 1) * N:
            tsamp = samples[(i * N) : (i + 1) * N]
            sampList.append(tsamp)
        else:
            if verbose:
                print("Fasta has only %d - not %d" % (len(samples), (i + 1) * N))

    return np.array(sampList)


def test_date_Ymd(x):
    x = str(x)
    try:
        _ = int(x[:4])
    except ValueError:
        return False

    try:
        _ = int(x[5:7])
    except ValueError:
        return False

    try:
        _ = int(x[8:])
    except ValueError:
        return False

    return True


def to_roman_numeral(value):

    roman_map = {
        1: "I",
        5: "V",
        10: "X",
        50: "L",
        100: "C",
        500: "D",
        1000: "M",
    }
    result = ""
    remainder = value

    for i in sorted(roman_map.keys(), reverse=True):
        if remainder > 0:
            multiplier = i
            roman_digit = roman_map[i]

            times = remainder // multiplier
            remainder = remainder % multiplier
            result += roman_digit * times

            if result.endswith("VIIII"):
                result = result.replace("VIIII", "IX")
            elif result.endswith("IIII"):
                result = result.replace("IIII", "IV")

    return result


def list_dfmeta_field(dfmeta, field):
    try:
        vals = dfmeta[field].unique()
    except ValueError:
        return []

    return np.sort(vals)


def create_dir(root: Path, subdir:str=''):

    if isinstance(root, str):
        root = Path(root)

    if isinstance(subdir, str) and subdir != "":
        path = root / subdir
    else:
        path = root

    os.makedirs(path, exist_ok=True)

    return path


def echo_print(stri: str, verbose: bool = False):
    if not isinstance(stri, str):
        return stri

    if verbose:
        print(stri)
    return stri + "\n"

def all_equal_list(cols1, cols2):
    cols1 = list(cols1)
    cols2 = list(cols2)

    if cols1 == [] and cols2 == []:
        return True

    if len(cols1) != len(cols2):
        return False

    soma = np.sum([1 if cols1[i] != cols2[i] else 0 for i in range(len(cols1))])
    return soma == 0


def zip_compress(zip_name, root_zip, root_data, type_of_file="pdf"):
    """
    https://stackoverflow.com/questions/47438424/python-zip-compress-multiple-files

    Select the compression mode ZIP_DEFLATED for compression
    or zipfile.ZIP_STORED to just store the file
    """
    compression = zipfile.ZIP_DEFLATED
    filefull = osjoin(root_zip, zip_name)
    print(f"Starting: {filefull}")

    files = [
        x
        for x in os.listdir(root_data)
        if x.endswith(f".{type_of_file}") and x != f".{type_of_file}"
    ]
    if len(files) == 0:
        print(f"There are no files in '{root_data}'")
        return False

    print(f"Saving '{filefull}' # {len(files)} ... {type_of_file} files...")
    # create the zip file first parameter path/name, second mode
    zf = zipfile.ZipFile(filefull, mode="w")

    try:
        ret = True
        i = -1
        for fname in files:
            i += 1
            if i % 1000 == 0:
                print(".", end="")
            zf.write(osjoin(root_data, fname), fname, compress_type=compression)
    except ValueError:
        print(f"\nError '{ValueError}', could not zip the files.")
        ret = False
    finally:
        # Don't forget to close the file!
        print(f"\nfile saved at '{filefull}' with {len(files)} files.")
        zf.close()

    return ret


def copy_file(
    fullname_ori: str, fullname_dest: str, force: bool = False, verbose: bool = False
) -> bool:

    if not exists(fullname_ori):
        print(f"File ORI does not exist: {fullname_ori}")
        return False
    if not force and exists(fullname_dest):
        if verbose:
            print(f"File DEST already exists: {fullname_dest}")
        return False

    try:
        shutil.copyfile(fullname_ori, fullname_dest)
    except ValueError:
        print(f"Error: {fullname_ori} -> {fullname_dest}")
        return False

    return True


def rename_file(root: str, file_ori: str, file_dest: str, verbose: bool = False) -> bool:

    filename_ori = osjoin(root, file_ori)
    filename_dest = osjoin(root, file_dest)

    if not exists(filename_ori):
        print(f"File ORI does not exist: {filename_ori}")
        return False
    if exists(filename_dest):
        if verbose:
            print(f"File DEST already exists: {filename_dest}")
        return False

    try:
        os.rename(filename_ori, filename_dest)
    except ValueError:
        print(f"Error: {filename_ori} -> {filename_dest}")
        return False

    return True


def rename_files(
    root, pattern_src: str, pattern_dst: str, _type: str = ".tsv", verbose: bool = False
) -> bool:

    if _type is None:
        lista = [x for x in os.listdir(root) if pattern_src in x and "~lock" not in x]
    else:
        lista = [
            x
            for x in os.listdir(root)
            if pattern_src in x and "~lock" not in x and x.endswith(_type)
        ]

    if lista == []:
        if verbose:
            print(f"There are no files in {root}")
        return False

    ret = 1
    for fname in lista:
        fname_src = osjoin(root, fname)

        fname_dst = fname.replace(pattern_src, pattern_dst)
        fname_dst = osjoin(root, fname_dst)

        if not exists(fname_src):
            print(f"fname does not exist: '{fname_src}'", fname_src)
            continue

        if exists(fname_dst):
            print(f"fname already exists: '{fname_dst}'")
            continue

        if fname_dst == fname_src:
            print("fname did not change: 'fname_src'")
            continue

        if verbose:
            print(fname_src, "to", fname_dst)

        try:
            os.rename(fname_src, fname_dst)
        except ValueError:
            ret *= 0

    return ret == 1


def series_round_scientific(df_series):
    mat_list = [x if isinstance(x, list) else eval(x) for x in df_series]

    for i in range(len(mat_list)):
        mat = mat_list[i]
        mat = [f"{float(x):.2e}" for x in mat]
        mat_list[i] = mat

    return mat_list


def which_semester(month: int) -> int:
    if not isint(month):
        return -1

    month = int(month)

    if month >= 1 and month <= 6:
        return 1

    if month >= 7 and month <= 12:
        return 2

    return -1


def df_to_md_table(df: pd.DataFrame, cols: List, rows: List, row_name: str, nround: int = 2) -> str:
    df = df.reset_index(drop=True).copy()

    head_row = f"| {row_name} | " + " | ".join(cols) + " |"
    sep_row = "| " + " | ".join(["---"] * (len(cols) + 1)) + " |"

    data_rows = []
    for i in range(len(df)):
        row_values = [str(np.round(val, nround)) for val in df.iloc[i][cols]]
        data_rows.append("| " + rows[i] + " | " + " | ".join(row_values) + " |")

    # Combine all parts
    markdown_table = "\n".join([head_row, sep_row] + data_rows)

    return markdown_table


def hex_to_rgb(value):
    value = value.lstrip("#")
    lv = len(value)
    return tuple(int(value[i : i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_hex(r: int, g: int, b: int) -> str:
    """
    Converts RGB color values to a hexadecimal string.

    Args:
            r (int): Red component (0-255).
            g (int): Green component (0-255).
            b (int): Blue component (0-255).

    Returns:
            str: Hexadecimal color string (e.g., "#RRGGBB").
    """
    return f"#{r:02x}{g:02x}{b:02x}"


def inc_rgb_to_hex(r: int, g: int, b: int, delta: int = 5, verbose: bool = False) -> str:
    """
    Converts RGB color values to a hexadecimal string.

    Args:
            r (int): Red component (0-255).
            g (int): Green component (0-255).
            b (int): Blue component (0-255).

    Returns:
            str: Hexadecimal color string (e.g., "#RRGGBB").
    """

    if r + delta <= 255:
        r += delta
    elif g + delta <= 255:
        g += delta
    elif b + delta <= 255:
        b += delta
    else:
        if verbose:
            print("RGB problems to increment")
        pass

    hex_color = f"#{r:02x}{g:02x}{b:02x}"
    # print("###", hex_color, r, g, b, delta)

    return hex_color


def text_starts_with_word(text: str, word: str) -> bool:

    pattern = r"^" + word
    ret = False
    for line in text.split("\n"):
        if re.match(pattern, line):
            ret = True
            break
    return ret


def create_empty_df(df3: pd.DataFrame, vals: list = []) -> pd.DataFrame:
    cols = list(df3.columns)

    if vals == []:
        vals = ["....."]
        vals += [None] * (len(cols) - 1)
    row = pd.Series(vals, index=cols)
    df_empty = pd.DataFrame([row.tolist()], columns=row.index)

    return df_empty
