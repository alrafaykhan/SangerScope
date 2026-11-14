"""Variant records (SNVs and indels) with HGVS naming and protein-level consequences."""

from __future__ import annotations

from functools import cmp_to_key
from typing import TYPE_CHECKING, Iterable, Iterator

from .symbols import GAP, UntranslatableError, aa_from_triple, jget, jsub

if TYPE_CHECKING:  # pragma: no cover
    from .analysis import Analysis

DELETION = -1
DELINS = 0
INSERTION = 1
DUPLICATION = 2

TYPE_NAMES = {DELETION: "deletion", DELINS: "delins", INSERTION: "insertion", DUPLICATION: "duplication"}


def _cmp(a, b) -> int:
    return (a > b) - (a < b)


class Variant:
    """Common fields of SNVs and indels. Indices are 1-based."""

    kind = "variant"

    def __init__(self) -> None:
        self.direction = 1
        self.c_index = ""
        self.hgvs = ""
        self.aa_change = ""
        self.alignment_index = 0
        self.fwd_trace_index = 0
        self.rev_trace_index = 0
        self.fwd_trace_char = ""
        self.rev_trace_char = ""
        self.hit_count = 1
        self.g_index = 0
        self.on_target = False
        self.zygosity = ""
        self.combined_expression = ""
        # filled in by sangerscope.quality.assess
        self.quality = "good"
        self.quality_notes: list[str] = []
        self.fwd_phred: int | None = None
        self.rev_phred: int | None = None

    # ---------------------------------------------------------------- ordering

    def compare(self, other: "Variant") -> int:
        """Display order; 0 means the two calls describe the same variant."""
        if isinstance(self, Indel):
            if isinstance(other, Indel):
                if self.zygosity != other.zygosity:
                    return _cmp(self.zygosity, other.zygosity)
                if self.g_index != other.g_index:
                    return self.g_index - other.g_index
                if self.hgvs != other.hgvs:
                    return _cmp(self.hgvs, other.hgvs)
            else:
                return -1
        elif isinstance(other, Indel):
            return 1

        if self.on_target == other.on_target:
            if self.g_index != other.g_index:
                return self.g_index - other.g_index
            if self.hgvs != other.hgvs:
                return _cmp(self.hgvs, other.hgvs)
            return _cmp(self.zygosity, other.zygosity)
        return -1 if self.on_target else 1

    # ----------------------------------------------------------------- display

    @property
    def protein(self) -> str:
        return self.aa_change if self.aa_change else "p.?"

    @property
    def description(self) -> str:
        return f"{self.hgvs}, {self.protein}"

    @property
    def detected_in(self) -> str:
        if self.hit_count == 2:
            return "Fwd, Rev"
        return "Fwd" if self.direction == 1 else "Rev"

    @property
    def equivalent_expressions(self) -> str:
        return self.combined_expression

    def to_dict(self) -> dict:
        return {
            "variant": self.hgvs,
            "protein": self.protein,
            "zygosity": self.zygosity,
            "type": self.kind,
            "detected_in": self.detected_in,
            "hit_count": self.hit_count,
            "on_target": self.on_target,
            "reference_position": self.g_index,
            "equivalent_expressions": self.equivalent_expressions.strip().rstrip(";").rstrip(", "),
            "quality": self.quality,
            "quality_notes": "; ".join(self.quality_notes),
            "forward_phred": "" if self.fwd_phred is None else self.fwd_phred,
            "reverse_phred": "" if self.rev_phred is None else self.rev_phred,
        }


class SNV(Variant):
    kind = "SNV"

    def __init__(
        self,
        ctx: "Analysis",
        ref_char: str,
        fwd_trace_char: str,
        rev_trace_char: str,
        direction: int,
        c_index: str,
        alignment_index: int,
        fwd_trace_index: int,
        rev_trace_index: int,
        coding: bool,
        g_index: int,
        on_target: bool,
        zygosity: str,
    ):
        super().__init__()
        self.ref_char = ref_char
        self.fwd_trace_char = fwd_trace_char
        self.rev_trace_char = rev_trace_char
        self.direction = direction
        self.c_index = c_index
        self.alignment_index = alignment_index
        self.fwd_trace_index = fwd_trace_index
        self.rev_trace_index = rev_trace_index
        self.g_index = g_index
        self.on_target = on_target
        self.zygosity = zygosity
        self.hgvs = f"{c_index}{ref_char}>{self.alt_char}"
        if coding:
            self._make_aa_change(ctx)

    @property
    def alt_char(self) -> str:
        return self.fwd_trace_char if self.direction == 1 else self.rev_trace_char

    def _make_aa_change(self, ctx: "Analysis") -> None:
        digits = "".join(ch for ch in self.c_index if ch.isascii() and ch.isdigit())
        if not digits:
            return
        c_number = int(digits)
        frame = (c_number - 1) % 3
        # The codon comes from the spliced coding sequence: a codon can be interrupted by an intron.
        original = ctx.reference.codon_at(self.g_index, frame)
        if original is None:  # not coding (e.g. 3' UTR) or codon not fully covered by the reference
            return
        mutated = original[:frame] + self.alt_char + original[frame + 1 :]
        codon_number = (c_number - 1) // 3 + 1
        aa = ctx.reference.aa
        try:
            self.aa_change = f"p.({aa(original)}{codon_number}{aa(mutated)})"
        except UntranslatableError:
            self.aa_change = "(untranslatable)"


class EquivExpression:
    """One of several equivalent descriptions of an indel in a repeated region."""

    __slots__ = ("g_index1", "g_index2", "hgvs")

    def __init__(self, g_index1: int, g_index2: int, hgvs: str):
        self.g_index1 = g_index1
        self.g_index2 = g_index2
        self.hgvs = hgvs


class Indel(Variant):
    kind = "indel"

    def __init__(
        self,
        ctx: "Analysis",
        zygosity: str,
        direction: int,
        indel_type: int,
        indel_start_index: int,
        indel_end_index: int,
        focused_index: int,
        indel_seq: str,
        on_target: bool,
    ):
        super().__init__()
        self._aps = ctx.aligned_points
        aps = self._aps
        self.direction = direction
        self.zygosity = zygosity
        self.type = indel_type
        self.indel_seq = indel_seq
        self.on_target = True
        # keyed by g_index1; kept in descending order of g_index1 when listed
        self._equiv: dict[int, EquivExpression] = {}

        focus = jget(aps, focused_index - 1)
        self.alignment_index = focused_index
        self.fwd_trace_index = focus.fwd_trace_index
        self.rev_trace_index = focus.rev_trace_index
        if direction == 1:
            self.fwd_trace_char = focus.fwd_char
            self.rev_trace_char = GAP
        elif direction == -1:
            self.fwd_trace_char = GAP
            self.rev_trace_char = focus.rev_char

        ap1 = jget(aps, indel_start_index - 1)
        ap2 = jget(aps, indel_end_index - 1)
        self.c_index = ap1.c_index
        self.c_index2 = ap2.c_index
        self.coding1 = ap1.coding
        self.coding2 = ap2.coding
        self.g_index = ap1.g_index
        self.g_index2 = ap2.g_index

        self._add_equiv(self._make_hgvs(indel_start_index, indel_end_index))

        # Shift the indel left and right while the resulting sequence stays the same.
        try:
            original = self._mutated_seq(self.g_index, self.g_index2)
            for step in (-1, 1):
                offset = step
                while True:
                    mutated = self._mutated_seq(self.g_index + offset, self.g_index2 + offset)
                    if mutated is None:
                        break
                    if mutated == original:
                        start = self._aligned_index_from_g_index(self.g_index + offset)
                        end = self._aligned_index_from_g_index(self.g_index2 + offset)
                        self._add_equiv(self._make_hgvs(start, end))
                    offset += step
        except IndexError:
            pass

        if self.coding1 and self.coding2 and ctx.formatter.first_number == 1:
            self._make_aa_change(ctx)

        # HGVS recommends the most 3' position
        right_most = None
        right_most_index = -1
        for expression in self.equiv_expressions:
            if right_most_index < expression.g_index2:
                right_most_index = expression.g_index2
                right_most = expression
        self.hgvs = right_most.hgvs

    # ---------------------------------------------------- equivalent expressions

    def _add_equiv(self, expression: EquivExpression) -> None:
        self._equiv.setdefault(expression.g_index1, expression)

    @property
    def equiv_expressions(self) -> list[EquivExpression]:
        """Equivalent expressions, most 3' first."""
        return [self._equiv[k] for k in sorted(self._equiv, reverse=True)]

    @property
    def equivalent_expressions(self) -> str:
        expressions = self.equiv_expressions
        out = []
        for i, expression in enumerate(expressions):
            if expression.hgvs == self.hgvs:
                continue
            out.append(expression.hgvs)
            if i != len(expressions) - 1:
                out.append("; ")
        text = "".join(out)
        if len(self.combined_expression) >= 1:
            text = self.combined_expression + "; " + text
        return text

    def _aligned_index_from_g_index(self, g_index: int) -> int:
        for i, ap in enumerate(self._aps):
            if ap.g_index == g_index:
                return i + 1
        return 0

    def _mutated_seq(self, g1: int, g2: int) -> str | None:
        """Reference of the aligned region after applying this indel at g1..g2."""
        aps = self._aps
        if g1 < aps[0].g_index or g2 > aps[-1].g_index:
            return None
        out = []
        dup = []
        for ap in aps:
            if ap.ref_char == GAP:
                continue
            g = ap.g_index
            if self.type == DELETION:
                if not (g1 <= g <= g2):
                    out.append(ap.ref_char)
            elif self.type == INSERTION:
                out.append(ap.ref_char)
                if g == g1:
                    out.append(self.indel_seq)
            elif self.type == DUPLICATION:
                if g1 <= g <= g2:
                    dup.append(ap.ref_char)
                out.append(ap.ref_char)
                if g == g2:
                    out.append("".join(dup))
            elif self.type == DELINS:
                if not (g1 <= g <= g2):
                    out.append(ap.ref_char)
                if g == g1:
                    out.append(self.indel_seq)
        return "".join(out)

    def _make_hgvs(self, aligned_index1: int, aligned_index2: int) -> EquivExpression:
        aps = self._aps
        ap1 = jget(aps, aligned_index1 - 1)
        ap2 = jget(aps, aligned_index2 - 1)
        c1 = ap1.c_index
        c2 = ap2.c_index
        c2_short = c2[2:] if len(c2) > 2 else c2

        hgvs = ""
        if self.type == DELETION:
            hgvs = f"{c1}del" if c1 == c2 else f"{c1}_{c2_short}del"
        elif self.type == INSERTION:
            hgvs = f"{c1}_{c2_short}ins{self.indel_seq}"
        elif self.type == DUPLICATION:
            hgvs = f"{c1}dup" if c1 == c2 else f"{c1}_{c2_short}dup"
        elif self.type == DELINS:
            hgvs = f"{c1}delins{self.indel_seq}" if c1 == c2 else f"{c1}_{c2_short}delins{self.indel_seq}"
        return EquivExpression(ap1.g_index, ap2.g_index, hgvs)

    # ------------------------------------------------------------ protein level

    def _make_aa_change(self, ctx: "Analysis") -> None:
        reference = ctx.reference
        ref = reference.sequence
        original = reference.spliced() + ref[reference.cdna_end[-1] :]
        self._aa = reference.aa
        try:
            self._aa_change_from(original)
        except UntranslatableError:
            self.aa_change = "(untranslatable)"
        except (IndexError, ValueError):
            return

    def _aa_change_from(self, original: str) -> None:
        c1 = int("".join(ch for ch in self.c_index if ch.isascii() and ch.isdigit()))
        c2 = int("".join(ch for ch in self.c_index2 if ch.isascii() and ch.isdigit()))
        length = 0
        shifted = ""
        if self.type == DELETION:
            length = self.g_index2 - self.g_index + 1
            shifted = jsub(original, 0, c1 - 1) + jsub(original, c2)
        elif self.type == INSERTION:
            length = len(self.indel_seq)
            shifted = jsub(original, 0, c1) + self.indel_seq + jsub(original, c2 - 1)
        elif self.type == DUPLICATION:
            length = self.g_index2 - self.g_index + 1
            shifted = jsub(original, 0, c1 - 1) + jsub(original, c1 - 1, c2) + jsub(original, c1 - 1)
        elif self.type == DELINS:
            length = len(self.indel_seq) - (self.g_index2 - self.g_index + 1)
            shifted = jsub(original, 0, c1 - 1) + self.indel_seq + jsub(original, c2)

        if length % 3 != 0:
            self._frameshift(original, shifted)
        else:
            self._in_frame(original, shifted)

    def _frameshift(self, original: str, shifted: str) -> None:
        prefix = ""
        i = 0
        while i < len(shifted):
            if i + 3 > len(original):
                break
            original_aa = self._aa(jsub(original, i, i + 3))
            shifted_aa = self._aa(jsub(shifted, i, i + 3))
            # a stop codon as the first changed amino acid is a nonsense change, not a frameshift
            if shifted_aa == "*":
                self.aa_change = f"p.({original_aa}{i // 3 + 1}{shifted_aa})"
                return
            if original_aa != shifted_aa:
                prefix = f"p.({original_aa}{i // 3 + 1}{shifted_aa}fs"
                break
            i += 3

        terminal_found = False
        count = 0
        try:
            while i < len(shifted):
                if i + 3 > len(original):
                    break
                count += 1
                if self._aa(jsub(shifted, i, i + 3)) == "*":
                    terminal_found = True
                    break
                i += 3
        except IndexError:
            pass
        self.aa_change = f"{prefix}*{count})" if terminal_found else f"{prefix})"

    def _in_frame(self, original: str, shifted: str) -> None:
        original_aa = [self._aa(original[i : i + 3]) for i in range(0, len(original) - 2, 3)]
        shifted_aa = [self._aa(shifted[i : i + 3]) for i in range(0, len(shifted) - 2, 3)]
        if original_aa == shifted_aa:
            self.aa_change = "no amino acid change"
            return

        left = 1
        while left <= len(original_aa) and left <= len(shifted_aa):
            if original_aa[left - 1] != shifted_aa[left - 1]:
                break
            left += 1
        left -= 1  # 1-based position of the last amino acid shared from the N-terminus
        original_right = len(original_aa)
        shifted_right = len(shifted_aa)
        while original_right > left and shifted_right > left:
            if original_aa[original_right - 1] != shifted_aa[shifted_right - 1]:
                break
            original_right -= 1
            shifted_right -= 1

        if left == shifted_right:  # deletion
            text = f"p.({jget(original_aa, left)}{left + 1}"
            if left + 1 != original_right:
                text += f"_{jget(original_aa, original_right - 1)}{original_right}"
            text += "del)"
        elif left == original_right:  # insertion or duplication
            size = shifted_right - left
            dup = True
            for i in range(left + 1, shifted_right + 1):
                if jget(shifted_aa, i - 1) != jget(shifted_aa, i - 1 - size):
                    dup = False
                    break
            if dup:
                text = f"p.({jget(original_aa, left - 1 - (size - 1))}{left - (size - 1)}"
                if size > 1:
                    text += f"_{jget(original_aa, left - 1)}{left}"
                text += "dup)"
            else:
                text = f"p.({jget(original_aa, left - 1)}{left}"
                if left < len(original_aa):
                    text += f"_{jget(original_aa, left)}{left + 1}"
                text += "ins" + "".join(jget(shifted_aa, i - 1) for i in range(left + 1, shifted_right + 1)) + ")"
        else:  # delins
            text = f"p.({jget(original_aa, left)}{left + 1}"
            if left + 1 != original_right:
                text += f"_{jget(original_aa, original_right - 1)}{original_right}"
            text += "delins" + "".join(jget(shifted_aa, i - 1) for i in range(left + 1, shifted_right + 1)) + ")"
        self.aa_change = text

    def to_dict(self) -> dict:
        d = super().to_dict()
        d["type"] = TYPE_NAMES[self.type]
        return d


class VariantSet:
    """Ordered set of variants; calls that compare equal are stored once."""

    def __init__(self, variants: Iterable[Variant] = ()):
        self._items: list[Variant] = []
        for v in variants:
            self.add(v)

    def _find(self, variant: Variant) -> int:
        for i, existing in enumerate(self._items):
            if existing.compare(variant) == 0:
                return i
        return -1

    def __contains__(self, variant: Variant) -> bool:
        return self._find(variant) >= 0

    def add(self, variant: Variant) -> None:
        if self._find(variant) < 0:
            self._items.append(variant)
            self._items.sort(key=cmp_to_key(lambda a, b: a.compare(b)))

    def remove(self, variant: Variant) -> None:
        i = self._find(variant)
        if i >= 0:
            del self._items[i]

    def copy(self) -> "VariantSet":
        new = VariantSet()
        new._items = list(self._items)
        return new

    def __iter__(self) -> Iterator[Variant]:
        return iter(list(self._items))

    def __len__(self) -> int:
        return len(self._items)

    def __getitem__(self, index: int) -> Variant:
        return self._items[index]
