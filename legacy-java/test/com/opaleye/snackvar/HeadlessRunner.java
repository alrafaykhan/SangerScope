package com.opaleye.snackvar;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.TreeSet;
import java.util.Vector;

import com.opaleye.snackvar.mmalignment.AlignedPair;
import com.opaleye.snackvar.mmalignment.MMAlignment;
import com.opaleye.snackvar.reference.Reference;
import com.opaleye.snackvar.tools.SymbolTools;
import com.opaleye.snackvar.variants.Variant;
import com.opaleye.snackvar.variants.VariantCallerFilter;

/**
 * Runs the analysis without the user interface, for comparing results with the Python implementation.
 * Input : a manifest file with one case per line (id, reference FASTA, forward AB1 or "-", reverse AB1 or "-",
 * separated by tabs). Output : lines starting with "@@" on standard output.
 * The steps are those of RootController.handleRun(), with the suggested trimming applied without confirmation.
 */
public class HeadlessRunner {
	private static PrintStream out;

	private static GanseqTrace load(File file, RootController rc, boolean complement) throws Exception {
		GanseqTrace trace = new GanseqTrace(file, rc);
		if(trace.getSequenceLength()<30) throw new Exception("too short");
		int start = trace.getFrontTrimPosition();
		int end = trace.getTailTrimPosition();
		if(start >= end) {
			start = 0;
			end = trace.traceLength * GanseqTrace.traceWidth-1;
		}
		trace.makeTrimmedTrace(start, end, complement);
		return trace;
	}

	private static int identity(AlignedPair ap) {
		int score = 0;
		for(int i=0;i<ap.getAlignedString1().length();i++)
			if(ap.getAlignedString1().charAt(i)==ap.getAlignedString2().charAt(i)) score++;
		return score;
	}

	private static void align(RootController rc) throws Exception {
		MMAlignment mma = new MMAlignment(rc.gapOpenPenalty);
		String ref = rc.reference.getRefString();
		AlignedPair fwdAp = null, cFwdAp = null, revAp = null, cRevAp = null;
		boolean fwdReversed = false, revReversed = false;
		if(rc.fwdLoaded) {
			fwdAp = mma.localAlignment(ref, rc.trimmedFwdTrace.getSequence());
			cFwdAp = mma.localAlignment(ref, SymbolTools.getComplementString(rc.trimmedFwdTrace.getSequence()));
			fwdReversed = identity(fwdAp) < identity(cFwdAp);
		}
		if(rc.revLoaded) {
			revAp = mma.localAlignment(ref, rc.trimmedRevTrace.getSequence());
			cRevAp = mma.localAlignment(ref, SymbolTools.getComplementString(rc.trimmedRevTrace.getSequence()));
			revReversed = identity(revAp) < identity(cRevAp);
		}
		if(!rc.revLoaded && fwdReversed) {
			rc.trimmedRevTrace = (GanseqTrace)rc.trimmedFwdTrace.clone();
			rc.trimmedRevTrace.makeComplement();
			rc.revLoaded = true;
			rc.trimmedFwdTrace = null;
			rc.fwdLoaded = false;
			fwdAp = null;
			revAp = cFwdAp;
		}
		else if(!rc.fwdLoaded && revReversed) {
			rc.trimmedFwdTrace = (GanseqTrace)rc.trimmedRevTrace.clone();
			rc.trimmedFwdTrace.makeComplement();
			rc.fwdLoaded = true;
			rc.trimmedRevTrace = null;
			rc.revLoaded = false;
			revAp = null;
			fwdAp = cRevAp;
		}
		else if(fwdReversed && revReversed) {
			GanseqTrace temp = (GanseqTrace)rc.trimmedRevTrace.clone();
			rc.trimmedRevTrace = (GanseqTrace)rc.trimmedFwdTrace.clone();
			rc.trimmedRevTrace.makeComplement();
			rc.trimmedFwdTrace = temp;
			rc.trimmedFwdTrace.makeComplement();
			fwdAp = cRevAp;
			revAp = cFwdAp;
		}
		if(rc.fwdLoaded && !rc.revLoaded)
			rc.alignedPoints = rc.formatter.format2(fwdAp, rc.reference, rc.trimmedFwdTrace, 1);
		else if(!rc.fwdLoaded && rc.revLoaded)
			rc.alignedPoints = rc.formatter.format2(revAp, rc.reference, rc.trimmedRevTrace, -1);
		else
			rc.alignedPoints = rc.formatter.format3(fwdAp, revAp, rc.reference, rc.trimmedFwdTrace, rc.trimmedRevTrace);
	}

	private static void setRange(RootController rc) {
		Vector<AlignedPoint> aps = rc.alignedPoints;
		rc.fwdTraceStart = 0; rc.fwdTraceEnd = 0; rc.revTraceStart = 0; rc.revTraceEnd = 0;
		boolean fwdFound = false, revFound = false;
		for(int i=0;i<aps.size();i++) {
			AlignedPoint ap = aps.get(i);
			if(!fwdFound && ap.getFwdChar() != Formatter.gapChar) {
				rc.fwdTraceStart = i+1;
				rc.trimmedFwdTrace.setAlignedRegionStart(ap.getFwdTraceIndex());
				fwdFound = true;
			}
			if(!revFound && ap.getRevChar() != Formatter.gapChar) {
				rc.revTraceStart = i+1;
				rc.trimmedRevTrace.setAlignedRegionStart(ap.getRevTraceIndex());
				revFound = true;
			}
			if(revFound && fwdFound) break;
		}
		fwdFound = false; revFound = false;
		for(int i=aps.size()-1;i>=0;i--) {
			AlignedPoint ap = aps.get(i);
			if(!fwdFound && ap.getFwdChar() != Formatter.gapChar) {
				rc.fwdTraceEnd = i+1;
				rc.trimmedFwdTrace.setAlignedRegionEnd(ap.getFwdTraceIndex());
				fwdFound = true;
			}
			if(!revFound && ap.getRevChar() != Formatter.gapChar) {
				rc.revTraceEnd = i+1;
				rc.trimmedRevTrace.setAlignedRegionEnd(ap.getRevTraceIndex());
				revFound = true;
			}
			if(revFound && fwdFound) break;
		}
		rc.startRange = Integer.max(rc.fwdTraceStart, rc.revTraceStart);
		rc.endRange = Integer.min(rc.fwdTraceEnd, rc.revTraceEnd);
	}

	private static Vector<Variant> detectHeteroIndel(RootController rc) {
		Vector<Variant> list = new Vector<Variant>();
		Variant fwdIndel = null, revIndel = null;
		if(rc.fwdLoaded) {
			fwdIndel = new HeteroTrace(rc.trimmedFwdTrace, rc).detectHeteroIndel();
			if(fwdIndel != null) list.add(fwdIndel);
		}
		if(rc.revLoaded) {
			revIndel = new HeteroTrace(rc.trimmedRevTrace, rc).detectHeteroIndel();
			if(revIndel != null) {
				if(fwdIndel != null && fwdIndel.getHGVS().equals(revIndel.getHGVS())) {
					revIndel.setHitCount(2);
					revIndel.setFwdTraceChar(fwdIndel.getFwdTraceChar());
					revIndel.setFwdTraceIndex(fwdIndel.getFwdTraceIndex());
				}
				list.add(revIndel);
			}
		}
		return list;
	}

	private static void runCase(String id, String refPath, String fwdPath, String revPath) {
		try {
			RootController rc = new RootController();
			rc.reference = new Reference(new File(refPath), Reference.FASTA);
			if(!fwdPath.equals("-")) {
				rc.trimmedFwdTrace = load(new File(fwdPath), rc, false);
				rc.fwdLoaded = true;
			}
			if(!revPath.equals("-")) {
				rc.trimmedRevTrace = load(new File(revPath), rc, true);
				rc.revLoaded = true;
			}
			TreeSet<Variant> variants = null;
			while(true) {
				rc.formatter = new Formatter(1);
				align(rc);
				setRange(rc);
				VariantCallerFilter vcf = new VariantCallerFilter(rc, detectHeteroIndel(rc));
				variants = vcf.getVariantList();
				if(vcf.misAlignment(variants) && rc.gapOpenPenalty < 1000) {
					if(rc.gapOpenPenalty == RootController.defaultGOP) rc.gapOpenPenalty = 200;
					else rc.gapOpenPenalty += 200;
					continue;
				}
				break;
			}
			StringBuffer ref = new StringBuffer(), fwd = new StringBuffer(), rev = new StringBuffer();
			for(AlignedPoint ap : rc.alignedPoints) {
				ref.append(ap.getRefChar()); fwd.append(ap.getFwdChar()); rev.append(ap.getRevChar());
			}
			out.println("@@GOP\t" + id + "\t" + rc.gapOpenPenalty);
			out.println("@@ALN\t" + id + "\t" + rc.alignedPoints.get(0).getGIndex() + "\t" + ref + "\t" + fwd + "\t" + rev);
			for(Variant v : variants) {
				out.println("@@VAR\t" + id + "\t" + v.getHGVS() + "\t" + v.getAAChange() + "\t" + v.getZygosity() + "\t"
						+ v.getHitCount() + "\t" + v.getFromProperty() + "\t" + v.getEquivalentExpressionsProperty());
			}
		}
		catch(Throwable ex) {
			out.println("@@ERR\t" + id + "\t" + ex.getClass().getSimpleName());
		}
		out.println("@@END\t" + id);
	}

	public static void main(String[] args) throws Exception {
		out = System.out;
		// the analysis classes print progress information; keep standard output for the results
		System.setOut(new PrintStream(new java.io.OutputStream() { public void write(int b) {} }));
		System.setErr(new PrintStream(new java.io.OutputStream() { public void write(int b) {} }));
		try (BufferedReader reader = new BufferedReader(new FileReader(args[0]))) {
			String line;
			while((line = reader.readLine()) != null) {
				String[] f = line.split("\t");
				if(f.length < 4) continue;
				runCase(f[0], f[1], f[2], f[3]);
			}
		}
		out.flush();
	}
}
