import re
from pathlib import Path
from typing import Callable, Iterator, List, Optional, Union

from rxn.utilities.regex import capturing

from .rdf_reaction import RdfReaction


class RdfParsingError(RuntimeError):
    """Exception for RDF parsing errors."""


class InvalidBlock(RdfParsingError):
    """Exception raised when a block of an RDF file cannot be processed."""

    def __init__(self, lines: List[str]):
        self.lines = lines
        super().__init__("Invalid block:\n" + "\n".join(lines))


class IncompleteReaction(RdfParsingError):
    """Exception for incomplete reaction in RDF file."""


BLOCK_TYPE_REGEX = re.compile(r"\$" + capturing(r"\w+"))
RIREG_REGEX = re.compile(r"\$RFMT \$[RM]IREG " + capturing(r"\d+"))
DTYPE_REGEX = re.compile(r"\$DTYPE " + capturing(".*"))
DATUM_REGEX = re.compile(r"\$DATUM " + capturing(".*"))


class ParsedReaction:
    """
    Reaction under construction during parsing of an RDF file.
    """

    def __init__(self) -> None:
        self.rireg: Optional[int] = None
        self.n_precursors: Optional[int] = None
        self.n_products: Optional[int] = None
        self.mols: List[str] = []
        self.dtypes: List[str] = []
        self.datums: List[str] = []

    def handle_line_block(self, lines: List[str]) -> None:
        block_type = self.block_type(lines)
        if block_type == "RFMT":
            self.handle_rfmt(lines)
        elif block_type == "RXN":
            self.handle_rxn(lines)
        elif block_type == "MOL":
            self.handle_mol(lines)
        elif block_type == "DTYPE":
            self.handle_dtype(lines)
        elif block_type == "DATUM":
            self.handle_datum(lines)
        else:
            raise ValueError(f"Invalid block type: {block_type}")

    def handle_rfmt(self, lines: List[str]) -> None:
        if len(lines) != 1:
            raise InvalidBlock(lines)
        match = RIREG_REGEX.match(lines[0])
        if match is None:
            raise InvalidBlock(lines)
        self.rireg = int(match.group(1))

    def handle_rxn(self, lines: List[str]) -> None:
        if len(lines) != 5:
            raise InvalidBlock(lines)

        self.n_precursors, self.n_products = [int(x) for x in lines[4].split()]

    def handle_mol(self, lines: List[str]) -> None:
        self.mols.append("\n".join(lines[1:]))

    def handle_dtype(self, lines: List[str]) -> None:
        if len(lines) != 1:
            raise InvalidBlock(lines)

        match = DTYPE_REGEX.match(lines[0])
        if match is None:
            raise InvalidBlock(lines)
        self.dtypes.append(match.group(1))

    def handle_datum(self, lines: List[str]) -> None:
        match = DATUM_REGEX.match(lines[0])
        if match is None:
            raise InvalidBlock(lines)
        lines[0] = match.group(1)

        # special case: leave out $MFMT, otherwise the property will not be
        # valid MolBlocks for parsing with RDKit.
        if lines[0] == "$MFMT":
            lines = lines[1:]

        self.datums.append("\n".join(lines))

    def block_type(self, lines: List[str]) -> str:
        """Get the type of a block: RXN, DTYPE, DATUM, etc."""
        match = BLOCK_TYPE_REGEX.match(lines[0])
        if match is None:
            raise InvalidBlock(lines)
        return match.group(1)

    def to_reaction(self) -> RdfReaction:
        if self.rireg is None or self.n_precursors is None or self.n_products is None:
            raise IncompleteReaction()

        if len(self.mols) != self.n_precursors + self.n_products:
            raise RdfParsingError()

        precursors = self.mols[: self.n_precursors]
        products = self.mols[-self.n_products :]

        if len(self.dtypes) != len(self.datums):
            raise RdfParsingError()

        meta = {key: value for key, value in zip(self.dtypes, self.datums)}

        return RdfReaction(
            reactants=precursors,
            reagents=[],
            products=products,
            meta=meta,
            reaction_index=self.rireg,
        )


class RdfParser:
    """
    Custom parser for RDF files.
    """

    def __init__(self, filename: Union[Path, str], encoding: str = "latin-1"):
        """
        Args:
            filename: path to the RDF file to read.
            encoding: file encoding. Defaults to latin-1 because Thieme has such
                an encoding for several files.
        """
        self.filename = filename
        self.encoding = encoding

    def __iter__(self) -> Iterator[RdfReaction]:
        yield from self.iter_reactions()

    def iter_reactions(self) -> Iterator[RdfReaction]:
        block_iterator = self.iter_blocks()
        # Consume line with RDFILE
        _ = next(block_iterator)
        # Consume line with DATM
        _ = next(block_iterator)

        current_reaction = None
        for lines in block_iterator:
            if lines[0].startswith("$RFMT"):
                # A new reaction started.
                # We yield the current one before initializing the new one
                if current_reaction is not None:
                    yield current_reaction.to_reaction()
                current_reaction = ParsedReaction()

            if current_reaction is None:
                raise RuntimeError("No reaction block started")

            current_reaction.handle_line_block(lines)

        # yield the last reaction
        if current_reaction is not None:
            yield current_reaction.to_reaction()

    def iter_blocks(self) -> Iterator[List[str]]:
        current_line_block: List[str] = []
        with open(self.filename, "rt", encoding=self.encoding) as f:
            for line in f:
                line = line.rstrip("\n")

                if line.startswith("$"):
                    if current_line_block:
                        yield current_line_block
                        current_line_block = []

                current_line_block.append(line)

        # Last line block at the end of the file
        if current_line_block:
            yield current_line_block


def iterate_reactions_from_file(
    filename: Union[Path, str],
    filter_fn: Optional[Callable[[RdfReaction], bool]] = None,
) -> Iterator[RdfReaction]:
    parser = RdfParser(filename)
    reactions = (entry for entry in parser.iter_reactions())

    if filter_fn is not None:
        reactions = (reaction for reaction in reactions if filter_fn(reaction))

    yield from reactions
