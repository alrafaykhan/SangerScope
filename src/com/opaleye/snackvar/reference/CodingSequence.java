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

package com.opaleye.snackvar.reference;

import java.util.List;

/**
 * Title : CodingSequence
 * Helpers for working on the spliced coding sequence of a genomic reference,
 * where the coding regions are given as 1-based, inclusive [start, end] intervals.
 */
public final class CodingSequence {

	private CodingSequence() {}

	/**
	 * Returns the 0-based offset of a genomic position within the spliced coding sequence.
	 * @param cDnaStart : 1-based start positions of the coding regions
	 * @param cDnaEnd : 1-based (inclusive) end positions of the coding regions
	 * @param gIndex : 1-based position on the reference
	 * @return offset in the spliced coding sequence, or -1 if gIndex is not in a coding region
	 */
	public static int splicedOffset(List<Integer> cDnaStart, List<Integer> cDnaEnd, int gIndex) {
		if(cDnaStart == null || cDnaEnd == null) return -1;

		int offset = 0;
		for(int i=0;i<cDnaStart.size();i++) {
			int start = cDnaStart.get(i).intValue();
			int end = cDnaEnd.get(i).intValue();
			if(gIndex < start) return -1;
			if(gIndex <= end) return offset + (gIndex - start);
			offset += (end - start + 1);
		}
		return -1;
	}

	/**
	 * Returns the coding regions of the reference joined together (introns and UTRs removed).
	 */
	public static String splice(String refString, List<Integer> cDnaStart, List<Integer> cDnaEnd) {
		StringBuilder sb = new StringBuilder();
		for(int i=0;i<cDnaStart.size();i++) {
			sb.append(refString, cDnaStart.get(i).intValue()-1, cDnaEnd.get(i).intValue());
		}
		return sb.toString();
	}

	/**
	 * Returns the codon containing a genomic position, taken from the spliced coding sequence
	 * so that codons interrupted by an intron are assembled from both flanking exons.
	 * @param refString : reference sequence
	 * @param cDnaStart : 1-based start positions of the coding regions
	 * @param cDnaEnd : 1-based (inclusive) end positions of the coding regions
	 * @param gIndex : 1-based position on the reference
	 * @param frame : position of gIndex within its codon (0, 1 or 2)
	 * @return the three bases of the codon,
	 *         or null if gIndex is not coding or the codon is not completely covered by the reference
	 */
	public static String codonAt(String refString, List<Integer> cDnaStart, List<Integer> cDnaEnd, int gIndex, int frame) {
		if(frame < 0 || frame > 2) return null;

		int offset = splicedOffset(cDnaStart, cDnaEnd, gIndex);
		if(offset < 0) return null;

		String spliced = splice(refString, cDnaStart, cDnaEnd);
		int codonStart = offset - frame;
		if(codonStart < 0 || codonStart + 3 > spliced.length()) return null;

		return spliced.substring(codonStart, codonStart + 3);
	}
}
