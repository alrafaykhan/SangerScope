/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.opaleye.snackvar.variants;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.RNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SymbolList;

import com.opaleye.snackvar.RootController;
import com.opaleye.snackvar.reference.CodingSequence;
import com.opaleye.snackvar.reference.Reference;

import javafx.beans.property.SimpleStringProperty;

/**
 * Title : SNV
 * A subclass of Variant
 * A class for SNV
 * @author Young-gon Kim
 *2018.10
 */
public class SNV extends Variant{
	private char refChar;

	/**
	 * 
	 * 	Constructor for a SNV
	 * @param refChar
	 * @param fwdTraceChar
	 * @param revTraceChar
	 * @param direction
	 * @param cIndex
	 * @param alignmentIndex
	 * @param fwdTraceIndex
	 * @param revTraceIndex
	 * @param coding
	 * @param gIndex
	 * @param hitCount
	 */
	
	public SNV(RootController rootController, char refChar, char fwdTraceChar, char revTraceChar, int direction, String cIndex, int alignmentIndex, int fwdTraceIndex, int revTraceIndex, boolean coding, int gIndex, boolean onTarget, String zygosity) {
		//tracePane adjust할때 fwdchar, revchar 모두 필요함. (gap이면 색칠안해야하니까)
		super();
		this.rootController = rootController;
		this.refChar = refChar;
		this.fwdTraceChar = fwdTraceChar;
		this.revTraceChar = revTraceChar;
		this.direction = direction;
		this.cIndex = cIndex;
		this.alignmentIndex = alignmentIndex;
		this.fwdTraceIndex = fwdTraceIndex;
		this.revTraceIndex = revTraceIndex;
		this.gIndex = gIndex;
		this.onTarget = onTarget;
		this.zygosity = zygosity;
		if(direction == 1) 
			HGVS = cIndex + refChar + ">" + fwdTraceChar;
		else
			HGVS = cIndex + refChar + ">" + revTraceChar;
		if(coding) makeAAChange();
		
		makeTableViewProperties();
		
	}

	/**
	 * A Getter for refChar 
	 */
	public char getRefChar() {
		return refChar;
	}

	/**
	 * Amino Acid change generation
	 */
	private void makeAAChange() {
		//gIndex : starts from 1
		Reference reference = rootController.reference;
		String originalSeq = null;
		String variantSeq = null;
		char variantChar = ' ';
		int i_cIndex = 0;
		int frameIndex = 0;
		if(direction == 1) 
			variantChar =  fwdTraceChar;
		else
			variantChar =  revTraceChar;
		try {
			i_cIndex = Integer.parseInt(cIndex.replaceAll("[^0-9]",""));
			// position of the mutated base within its codon (0, 1, 2)
			int frame = (i_cIndex-1) % 3;

			// The codon is taken from the spliced coding sequence, not from adjacent bases of the
			// genomic reference, because a codon can be interrupted by an intron.
			originalSeq = CodingSequence.codonAt(reference.getRefString(), reference.getcDnaStart(), reference.getcDnaEnd(), gIndex, frame);

			// not in a coding region (e.g. 3' UTR), or codon not completely covered by the reference
			if(originalSeq == null) return;

			variantSeq = originalSeq.substring(0, frame) + variantChar + originalSeq.substring(frame+1);
			frameIndex = (i_cIndex-1)/3+1;
		}
		catch (Exception ex) {
			ex.printStackTrace();
			return;
		}

		String originalAA = "";
		String variantAA = "";
		try {
			originalAA = Variant.getAAfromTriple(originalSeq);
			variantAA = Variant.getAAfromTriple(variantSeq);
			AAChange = "p.(" + originalAA + frameIndex + variantAA+ ")";
			//HGVS = HGVS + ", " + AAChange;
		}
		catch(IllegalAlphabetException iae) {
			//probably due to 'N' from indelSeq
			AAChange = "(untranslatable)";
			return;
			
		}
		catch (Exception ex) {
			ex.printStackTrace();
			return;
		}
		
	}
}
