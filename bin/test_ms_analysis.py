#!/usr/bin/env python

import sys
from pathlib import Path
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import ms_analysis as msa


def test_head3_path_only(h3_path: str):
    """If h3_path is not a file path for head3.lst,
    add file name and test for existence."""

    h3 = Path(h3_path)
    if h3.is_dir():
        h3 = h3.joinpath("head3.lst")

    if h3.is_file():
        if h3.name == "head3.lst":
            print(f"head3 file from path only exists: {h3}")
            return True
        else:
            print(f"Unexpected head3 file name: {h3.name}")
            return False
    else:
        print(f"test_head3_path_only - File not found: {h3}")
        return False
    
def test_alt_head3_path(h3_path):

    h3 = Path(h3_path)
    if not h3.is_file():
        print(f"Test expects a file; not found: {h3}.")
        return False

    pw = Path.cwd()
    if h3 == Path("head3.lst").resolve():
        print(f"Expected outcome using '{h3}':\n{len(msa.conformers) = }. Should be > 0: {len(msa.conformers)>0}")
        assert len(msa.conformers) > 0
    else:
        # redo conformer loading with non-standard path:
        new_conformers = msa.read_conformers(h3)
        n_new = len(new_conformers)
        
        if h3.parent == pw:
            print(f"Expected outcome using '{h3}':\n{len(msa.conformers) = }. Should be > 0: {len(msa.conformers)>0}")
            assert len(msa.conformers) > 0
        else:
            print(f"Expected outcome using len(msa.read_conformers({h3}) = {n_new}. Should be > 0: {n_new>0}")
            assert n_new > 0


def cli_parser():

    def resolved_path(p: str):
        """Resolved the command line path."""
        return Path(p).resolve()

    p = ArgumentParser(
        prog="test_ms_analysis",
        description="Provide unit tests on added features in bin/ms_analysis.py.",
        formatter_class=RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "-test_alt_head3_path",
        type=resolved_path,
        default=None,
        help="""Test the difference of the conformers list loaded from
        ms_analysis.py (using default head3 path) vs that created when an alternate
        path to head3.lst is provided at the command line.""",
    )

    return p


def main(argv=None):
    """Command line interface to test added features in bin/ms_analysis.py."""

    cli_parse = cli_parser()
    args = cli_parse.parse_args(argv)

    if args.test_alt_head3_path is None or argv is None:
        cli_parse.print_help()
        return
    
    h3 = Path(args.test_alt_head3_path)
    print(f"Resolved cli path: {h3}\npwd: {Path.cwd()}")

    if h3.name == "head3.lst":
        test_head3_path_only(h3.parent)
    elif h3.is_dir():
        test_head3_path_only(h3)
    
    if h3.is_file():
        test_alt_head3_path(h3)


if __name__ == "__main__":
    #sys.exit(main(sys.argv[1:]))
    main(sys.argv[1:])
