package com.opaleye.snackvar.reference;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;

import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.Test;

class CodingSequenceTest {

	// 5'UTR  exon1   intron   exon2    3'UTR
	// gg     ATGGC   gtaag    TTGGTAA  cc
	// 1-2    3-7     8-12     13-19    20-21
	// CDS : ATG GCT TGG TAA  -> codon 2 (GCT) is split as GC | T by the intron
	private static final String REF = "ggATGGCgtaagTTGGTAAcc";
	private static final List<Integer> START = Arrays.asList(3, 13);
	private static final List<Integer> END = Arrays.asList(7, 19);

	@Test
	void splicesCodingRegions() {
		assertEquals("ATGGCTTGGTAA", CodingSequence.splice(REF, START, END));
	}

	@Test
	void splicedOffsetSkipsIntrons() {
		assertEquals(0, CodingSequence.splicedOffset(START, END, 3));
		assertEquals(4, CodingSequence.splicedOffset(START, END, 7));
		assertEquals(5, CodingSequence.splicedOffset(START, END, 13));
		assertEquals(11, CodingSequence.splicedOffset(START, END, 19));
	}

	@Test
	void nonCodingPositionsHaveNoOffset() {
		assertEquals(-1, CodingSequence.splicedOffset(START, END, 2));	// 5' UTR
		assertEquals(-1, CodingSequence.splicedOffset(START, END, 10));	// intron
		assertEquals(-1, CodingSequence.splicedOffset(START, END, 20));	// 3' UTR
	}

	@Test
	void codonWithinOneExon() {
		assertEquals("ATG", CodingSequence.codonAt(REF, START, END, 3, 0));
		assertEquals("ATG", CodingSequence.codonAt(REF, START, END, 5, 2));
		assertEquals("TGG", CodingSequence.codonAt(REF, START, END, 15, 1));
	}

	@Test
	void codonSplitByIntronIsAssembledFromBothExons() {
		// c.4 and c.5 are the last two bases of exon 1, c.6 is the first base of exon 2
		assertEquals("GCT", CodingSequence.codonAt(REF, START, END, 6, 0));
		assertEquals("GCT", CodingSequence.codonAt(REF, START, END, 7, 1));
		assertEquals("GCT", CodingSequence.codonAt(REF, START, END, 13, 2));
	}

	@Test
	void noCodonOutsideCodingRegion() {
		assertNull(CodingSequence.codonAt(REF, START, END, 10, 0));	// intron
		assertNull(CodingSequence.codonAt(REF, START, END, 20, 0));	// 3' UTR
	}

	@Test
	void noCodonWhenReferenceDoesNotCoverIt() {
		// partial reference whose coding region starts in the middle of a codon
		assertNull(CodingSequence.codonAt("GCTTGG", Arrays.asList(1), Arrays.asList(6), 1, 2));
		assertNull(CodingSequence.codonAt("ATGGC", Arrays.asList(1), Arrays.asList(5), 4, 0));
	}
}
