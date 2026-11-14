package com.opaleye.snackvar.variants;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import com.opaleye.snackvar.RootController;
import com.opaleye.snackvar.reference.Reference;

class SNVTest {

	// 5'UTR  exon1   intron   exon2    3'UTR
	// gg     ATGGC   gtaag    TTGGTAA  cc
	// 1-2    3-7     8-12     13-19    20-21
	// CDS : ATG GCT TGG TAA  -> codon 2 (GCT, Ala) is split as GC | T by the intron
	private static final String FASTA = ">test\nggATGGCgtaagTTGGTAAcc\n";

	@TempDir
	Path tempDir;

	private RootController rootController;

	@BeforeEach
	void loadReference() throws Exception {
		File refFile = tempDir.resolve("test.fasta").toFile();
		Files.write(refFile.toPath(), FASTA.getBytes(StandardCharsets.US_ASCII));
		rootController = new RootController();
		rootController.reference = new Reference(refFile, Reference.FASTA);
	}

	private SNV forwardSnv(char refChar, char altChar, String cIndex, int gIndex) {
		return new SNV(rootController, refChar, altChar, ' ', 1, cIndex, 1, 1, 1, true, gIndex, true, "homo");
	}

	@Test
	void codonWithinOneExon() {
		// c.2T>C : ATG -> ACG
		assertEquals("p.(Met1Thr)", forwardSnv('T', 'C', "c.2", 4).getAAChange());
	}

	@Test
	void codonSplitByIntronBeforeTheJunction() {
		// c.5C>A, last base of exon 1 : GCT -> GAT.
		// Reading adjacent genomic bases would give GCg -> GAg (Glu).
		assertEquals("p.(Ala2Asp)", forwardSnv('C', 'A', "c.5", 7).getAAChange());
	}

	@Test
	void codonSplitByIntronAfterTheJunction() {
		// c.6T>G, first base of exon 2 : GCT -> GCG.
		// Reading adjacent genomic bases would give aaT -> aaG (Asn -> Lys).
		assertEquals("p.(Ala2Ala)", forwardSnv('T', 'G', "c.6", 13).getAAChange());
	}

	@Test
	void threePrimeUtrHasNoAminoAcidChange() {
		// shown as "p.?" rather than a change computed from non-coding bases
		assertEquals("c.*1C>T, p.?", forwardSnv('C', 'T', "c.*1", 20).getVariantProperty());
	}
}
