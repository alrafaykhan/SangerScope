"""Variant calling from the merged alignment, and the filters applied to the calls."""

from __future__ import annotations

from typing import TYPE_CHECKING

from .symbols import GAP, is_acgt, iupac_to_symbols, jget
from .variants import DELETION, DELINS, DUPLICATION, INSERTION, Indel, SNV, Variant, VariantSet

if TYPE_CHECKING:  # pragma: no cover
    from .analysis import Analysis

#: A single-strand call is dropped when the other strand shows the reference base at this quality or better.
FILTER_QUALITY_CUTOFF = 25


class VariantCaller:
    def __init__(self, ctx: "Analysis", hetero_indels: list[Indel]):
        self.ctx = ctx
        self.hetero_indels = list(hetero_indels)
        self.fwd_loaded = ctx.fwd_loaded
        self.rev_loaded = ctx.rev_loaded
        self.start_range = ctx.start_range
        self.end_range = ctx.end_range
        self.fwd = ctx.trimmed_fwd
        self.rev = ctx.trimmed_rev
        self.aps = ctx.aligned_points

    def call(self) -> VariantSet:
        variants = self.make_variant_list()
        variants = self._comparison_filter(variants)
        variants = self._compressed_peak_filter(variants)
        try:
            variants = self._make_homo_delins(variants)
        except IndexError:
            pass
        return variants

    def _on_target(self, index: int) -> bool:
        if not self.fwd_loaded or not self.rev_loaded:
            return True
        return self.start_range <= index <= self.end_range

    # ------------------------------------------------------------------ calling

    def make_variant_list(self) -> VariantSet:
        ctx = self.ctx
        aps = self.aps
        variants = VariantSet()

        # SNVs and homozygous indels are only looked for outside the double-peak region
        # caused by a heterozygous indel (0-based column where that region starts).
        fwd_limit = len(aps)
        rev_limit = -1
        for indel in self.hetero_indels:
            if indel.direction == 1:
                fwd_limit = indel.alignment_index - 1
            elif indel.direction == -1:
                rev_limit = indel.alignment_index - 1

        if len(self.hetero_indels) == 2 and self.hetero_indels[0].hgvs == self.hetero_indels[1].hgvs:
            self.hetero_indels = [self.hetero_indels[1]]
        for indel in self.hetero_indels:
            variants.add(indel)

        def snvs(i: int, direction: int) -> None:
            ap = aps[i]
            trace_char = ap.fwd_char if direction == 1 else ap.rev_char
            if is_acgt(trace_char):
                alts = [trace_char] if trace_char != ap.ref_char else []
                zygosity = "homo"
            else:
                alts = [b for b in iupac_to_symbols(trace_char) if b != ap.ref_char]
                zygosity = "hetero"
            for alt in alts:
                fwd_char = alt if direction == 1 else ap.fwd_char
                rev_char = alt if direction == -1 else ap.rev_char
                variant = SNV(
                    ctx, ap.ref_char, fwd_char, rev_char, direction, ap.c_index, i + 1,
                    ap.fwd_trace_index, ap.rev_trace_index, ap.coding, ap.g_index, self._on_target(i + 1), zygosity,
                )
                if direction == -1 and variant in variants:
                    variants.remove(variant)
                    variant.hit_count = 2
                variants.add(variant)

        for i in range(fwd_limit):
            ap = aps[i]
            if self.fwd_loaded and ap.fwd_char != GAP and ap.ref_char != GAP:
                snvs(i, 1)
        for i in range(rev_limit + 1, len(aps)):
            ap = aps[i]
            if self.rev_loaded and ap.rev_char != GAP and ap.ref_char != GAP:
                snvs(i, -1)

        if self.fwd_loaded:
            self._homo_insertions(variants, 1, 0, fwd_limit)
        if self.rev_loaded:
            self._homo_insertions(variants, -1, rev_limit + 1, len(aps))
        if self.fwd_loaded:
            self._homo_deletions(variants, 1, 0, fwd_limit)
        if self.rev_loaded:
            self._homo_deletions(variants, -1, rev_limit + 1, len(aps))
        return variants

    def _homo_insertions(self, variants: VariantSet, direction: int, begin: int, end: int) -> None:
        aps = self.aps

        def trace_char(ap) -> str:
            return ap.fwd_char if direction == 1 else ap.rev_char

        i = begin
        while i < end:
            ap = aps[i]
            if ap.ref_char == GAP and trace_char(ap) != GAP:
                j = i
                inserted = []
                while j < end:
                    ap2 = aps[j]
                    if ap2.ref_char == GAP and trace_char(ap2) != GAP:
                        inserted.append(trace_char(ap2))
                    else:
                        break
                    j += 1
                inserted_seq = "".join(inserted)

                duplication = False
                dup_start = i - len(inserted_seq)
                if dup_start >= 0:
                    duplication = all(
                        trace_char(aps[dup_start + k]) == inserted_seq[k] for k in range(len(inserted_seq))
                    )
                if duplication:
                    variant = Indel(
                        self.ctx, "homo", direction, DUPLICATION, dup_start + 1, dup_start + len(inserted_seq),
                        dup_start + 1, inserted_seq, self._on_target(dup_start + 1),
                    )
                else:
                    variant = Indel(
                        self.ctx, "homo", direction, INSERTION, i + 1, j + 1, i + 1, inserted_seq,
                        self._on_target(i + 1),
                    )
                if direction == -1 and variant in variants:
                    variants.remove(variant)
                    variant.hit_count = 2
                    variant.fwd_trace_char = ap.fwd_char
                variants.add(variant)
                i = j
            i += 1

    def _homo_deletions(self, variants: VariantSet, direction: int, begin: int, end: int) -> None:
        aps = self.aps

        def trace_char(ap) -> str:
            return ap.fwd_char if direction == 1 else ap.rev_char

        i = begin
        while i < end:
            ap = aps[i]
            if trace_char(ap) == GAP and ap.ref_char != GAP:
                j = i
                while j < end:
                    if not (trace_char(aps[j]) == GAP and ap.ref_char != GAP):
                        break
                    j += 1
                # a stretch at either end covered by only one strand is not a deletion
                if i == 0 or j == len(aps):
                    i = j + 1
                    continue
                variant = Indel(self.ctx, "homo", direction, DELETION, i + 1, j, i + 1, "", self._on_target(i + 1))
                if direction == -1 and variant in variants:
                    variants.remove(variant)
                    variant.hit_count = 2
                    variant.fwd_trace_char = ap.fwd_char
                variants.add(variant)
                i = j
            i += 1

    # ------------------------------------------------------------------ filters

    def _comparison_filter(self, variants: VariantSet) -> VariantSet:
        """Drop single-strand calls contradicted by a good-quality reference base on the other strand."""
        if not (self.fwd_loaded and self.rev_loaded):
            return variants
        aps = self.aps
        kept = variants.copy()

        def opposite(ap, direction: int) -> tuple[int, str]:
            if direction == 1:
                return ap.rev_quality, ap.rev_char
            if direction == -1:
                return ap.fwd_quality, ap.fwd_char
            return 0, "N"

        def contradicted(ap, direction: int) -> bool:
            quality, char = opposite(ap, direction)
            return quality >= FILTER_QUALITY_CUTOFF and char == ap.ref_char

        for v in variants:
            if v.hit_count != 1 or v.alignment_index < self.start_range or v.alignment_index > self.end_range:
                continue
            if isinstance(v, Indel) and v.zygosity == "hetero":
                continue
            if isinstance(v, SNV):
                if contradicted(jget(aps, v.alignment_index - 1), v.direction):
                    kept.remove(v)
            elif isinstance(v, Indel):
                if v.type == DELETION:
                    remove = True
                    counter = 0
                    while v.alignment_index - 1 + counter < len(aps):
                        ap = jget(aps, v.alignment_index - 1 + counter)
                        if not contradicted(ap, v.direction):
                            remove = False
                            break
                        if ap.g_index == v.g_index2:
                            break
                        counter += 1
                    if remove:
                        kept.remove(v)
                elif v.type == INSERTION:
                    if len(v.indel_seq) != 1:
                        continue
                    try:
                        before = jget(aps, v.alignment_index - 2)
                        after = jget(aps, v.alignment_index)
                    except IndexError:
                        continue
                    if contradicted(before, v.direction) and contradicted(after, v.direction):
                        kept.remove(v)
        return kept

    def _compressed_peak_filter(self, variants: VariantSet) -> VariantSet:
        """Drop the SNV + 1-base deletion pair produced when two neighbouring peaks are compressed into one."""
        filtered: list[str] = []
        formatter = self.ctx.formatter
        for loaded, trace, coordinate_map, direction in (
            (self.fwd_loaded, self.fwd, formatter.fwd_coordinate_map, 1),
            (self.rev_loaded, self.rev, formatter.rev_coordinate_map, -1),
        ):
            if not loaded:
                continue
            for i in range(trace.sequence_length):
                try:
                    self._compressed_peak_at(trace, coordinate_map, direction, i, filtered)
                except (IndexError, KeyError):
                    continue
        return VariantSet(v for v in variants if v.hgvs not in filtered)

    def _compressed_peak_at(self, trace, coordinate_map, direction: int, i: int, filtered: list[str]) -> None:
        aps = self.aps
        bases = iupac_to_symbols(trace.sequence[i])
        if len(bases) != 2:
            return
        base1, base2 = bases
        index = coordinate_map[i + 1]
        point = jget(aps, index - 1)
        if base1 == point.ref_char:
            ref_base, next_base = base1, base2
        elif base2 == point.ref_char:
            ref_base, next_base = base2, base1
        else:
            return
        second = trace.channel(next_base)
        position = jget(trace.base_calls, i)
        before, here, after = jget(second, position - 1), jget(second, position), jget(second, position + 1)
        if before < here < after:
            step = 1
        elif before > here > after:
            step = -1
        else:
            return
        index2 = index + step
        while True:
            point2 = jget(aps, index2 - 1)
            if point2.ref_char != next_base:
                break
            trace_char = point2.fwd_char if direction == 1 else point2.rev_char
            if trace_char == GAP:
                filtered.append(f"{point.c_index}{ref_base}>{next_base}")
                filtered.append(f"{point2.c_index}del{next_base}")
            index2 += step

    # ------------------------------------------------- combined delins expression

    def misalignment(self, variants: VariantSet) -> bool:
        """A homozygous indel next to a heterozygous indel usually means the aligner split the latter."""
        for v1 in variants:
            if not (isinstance(v1, Indel) and v1.zygosity == "hetero"):
                continue
            het_min = min(e.g_index1 for e in v1.equiv_expressions)
            het_max = max(e.g_index2 for e in v1.equiv_expressions)
            for v2 in variants:
                if not (isinstance(v2, Indel) and v2.zygosity == "homo"):
                    continue
                homo_min = min(e.g_index1 for e in v2.equiv_expressions) - 50
                homo_max = max(e.g_index2 for e in v2.equiv_expressions) + 50
                # the margin around the heterozygous indel grows with every homozygous indel examined
                het_min -= 50
                het_max += 50
                if not (het_max < homo_min or homo_max < het_min):
                    return True
        return False

    def _make_homo_delins(self, variants: VariantSet) -> VariantSet:
        """Annotate adjacent variants with the single delins expression they amount to."""
        ret = variants.copy()
        group = VariantSet()
        array = list(variants)
        building = False
        connected = False
        skip = False

        # adjacent homozygous variants
        for i in range(1, len(array)):
            if skip:
                skip = False
                continue
            v1, v2 = array[i - 1], array[i]
            v2_multiple = False
            if v1.zygosity == "hetero" or v2.zygosity == "hetero":
                connected = False
            else:
                v1_right = 0
                v2_left = 0
                if isinstance(v1, SNV):
                    v1_right = v1.g_index
                elif isinstance(v1, Indel):
                    v1_right = v1.equiv_expressions[-1].g_index2
                    if v1.type == INSERTION:
                        v1_right -= 1
                if isinstance(v2, SNV):
                    v2_left = v2.g_index
                elif isinstance(v2, Indel):
                    expressions = v2.equiv_expressions
                    v2_multiple = len(expressions) > 1
                    v2_left = expressions[0].g_index1
                    if v2.type == INSERTION:
                        v2_left += 1
                connected = v2_left <= v1_right + 1

            if building:
                if connected:
                    group.add(v2)
                    if v2_multiple:
                        ret = self._update_variant_list(ret, group)
                        group = VariantSet()
                        building = False
                        skip = True
                else:
                    ret = self._update_variant_list(ret, group)
                    group = VariantSet()
                    building = False
            elif connected:
                group.add(v1)
                group.add(v2)
                if v2_multiple:
                    ret = self._update_variant_list(ret, group)
                    group = VariantSet()
                    building = False
                    skip = True
                else:
                    building = True

        # consecutive heterozygous SNVs
        for i in range(1, len(array)):
            if skip:
                skip = False
                continue
            v1, v2 = array[i - 1], array[i]
            if v1.zygosity == "homo" or v2.zygosity == "homo" or isinstance(v1, Indel) or isinstance(v2, Indel):
                connected = False
            else:
                connected = v2.g_index == v1.g_index + 1
            if building:
                if connected:
                    group.add(v2)
                else:
                    ret = self._update_variant_list(ret, group)
                    group = VariantSet()
                    building = False
            elif connected:
                group.add(v1)
                group.add(v2)
                building = True
        return ret

    def _c_index_from_g_index(self, g_index: int) -> str:
        for ap in self.aps:
            if ap.g_index == g_index:
                return ap.c_index
        return ""

    def _aligned_index_from_g_index(self, g_index: int) -> int:
        for i, ap in enumerate(self.aps):
            if ap.g_index == g_index:
                return i + 1
        return 0

    def _update_variant_list(self, original: VariantSet, targets: VariantSet) -> VariantSet:
        """Apply all `targets` to the reference and describe the net change as one delins."""
        ret = original.copy()
        ref_chars: list[str] = []
        ref_index: list[int] = []
        for ap in self.aps:
            if ap.ref_char != GAP:
                ref_chars.append(ap.ref_char)
                ref_index.append(ap.g_index)
        ref_string = "".join(ref_chars)
        delins_string = ref_string
        delins_index = list(ref_index)

        target_list = list(targets)
        for v in target_list:
            index1 = -1
            index2 = 0 if isinstance(v, SNV) else -1
            for j, g in enumerate(delins_index):
                if g == v.g_index:
                    index1 = j
                if isinstance(v, Indel) and g == v.g_index2:
                    index2 = j
            if index1 == -1 or index2 == -1:
                return ret

            if isinstance(v, SNV):
                delins_string = delins_string[:index1] + v.alt_char + delins_string[index1 + 1 :]
            elif v.type == DELETION:
                delins_string = delins_string[:index1] + delins_string[index2 + 1 :]
                for _ in range(index1, index2 + 1):
                    del delins_index[index1]
            elif v.type == INSERTION:
                delins_string = delins_string[: index1 + 1] + v.indel_seq + delins_string[index2:]
                for _ in v.indel_seq:
                    delins_index.insert(index1 + 1, 0)
            elif v.type == DUPLICATION:
                delins_string = delins_string[: index2 + 1] + v.indel_seq + delins_string[index2 + 1 :]
                for _ in v.indel_seq:
                    delins_index.insert(index2 + 1, 0)

        left = 0
        original_right = len(ref_string) - 1
        shifted_right = len(delins_string) - 1
        while left < len(ref_string) and left < len(delins_string):
            if ref_string[left] != delins_string[left]:
                break
            left += 1
        left -= 1
        while original_right > left and shifted_right > left:
            if ref_string[original_right] != delins_string[shifted_right]:
                break
            original_right -= 1
            shifted_right -= 1

        combined = ""
        if left != shifted_right and left != original_right:
            g1 = jget(ref_index, left + 1)
            g2 = jget(ref_index, original_right)
            c1 = self._c_index_from_g_index(g1)
            c2 = self._c_index_from_g_index(g2)
            combined = c1 if c1 == c2 else f"{c1}_{c2[2:] if len(c2) > 2 else c2}"
            inserted = delins_string[left + 1 : shifted_right + 1]
            combined += "delins" + inserted
            start = self._aligned_index_from_g_index(g1)
            indel = Indel(self.ctx, "homo", 1, DELINS, start, self._aligned_index_from_g_index(g2), start, inserted, True)
            combined += ", " + indel.aa_change

        for v in target_list:
            ret.remove(v)
            v.combined_expression = combined
            ret.add(v)
        return ret
