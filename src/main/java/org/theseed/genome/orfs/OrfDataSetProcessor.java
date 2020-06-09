package org.theseed.genome.orfs;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.text.TextStringBuilder;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.io.ILabeledOutputStream;
import org.theseed.locations.Location;
import org.theseed.locations.SequenceLocation;
import org.theseed.proteins.CodonSet;
import org.theseed.proteins.DnaTranslator;
import org.theseed.utils.BaseProcessor;

public abstract class OrfDataSetProcessor extends BaseProcessor {

    // FIELDS

    protected static Logger log = LoggerFactory.getLogger(StartTrainProcessor.class);
    /** number of true starts output */
    private int trueCount;
    /** number of false starts output */
    private int falseCount;
    /** number of orfs processed */
    private int orfCount;
    /** output stream */
    private ILabeledOutputStream outStream;

    // COMMAND-LINE OPTIONS

    /** number of positions to output to the left of the candidate start */
    @Option(name = "--left", metaVar = "100", usage = "number of positions to output to the left of the candidate base pair")
    private int numLeft;

    /** number of positions to output to the right of the candidate start */
    @Option(name = "--right", metaVar = "50", usage = "number of positions to output to the right of the candidate base pair")
    private int numRight;

    /**
     * Initialize this object.
     */
    protected void setupDefaults() {
        this.numLeft = 52;
        this.numRight = 20;
        this.trueCount = 0;
        this.falseCount = 0;
        this.orfCount = 0;
    }

    /**
     * Validate the command-line options.
     */
    protected void validateCommonParms() {
        // Insure the left and right indicators are positive.
        if (this.numLeft < 0)
            throw new IllegalArgumentException("Left positions must be non-negative.");
        if (this.numRight < 0)
            throw new IllegalArgumentException("Right positions must be non-negative.");
    }

    /**
     * Write the header line to the output.
     */
    protected void writeHeader() {
        TextStringBuilder lineBuffer = new TextStringBuilder(6 * this.numLeft + 14 + 5 * this.numRight);
        lineBuffer.append("name");
        for (int i = -this.numLeft; i <= this.numRight; i++) {
            lineBuffer.appendSeparator('\t');
            lineBuffer.append("p.%d", i);
        }
        this.outStream.writeImmediate("type", lineBuffer.toString());
    }

    /**
     * Output start data for a single feature.
     *
     * @param genome	genome containing the feature
     * @param peg		protein feature whose start instances are to be output
     */
    protected void output(Genome genome, Feature peg) {
        // Get the location of the peg.
        Location pegLoc = peg.getLocation();
        // Extend it to an ORF.
        SequenceLocation orfLoc = pegLoc.createORF(genome);
        // Get the position of the original location's begin point in the sequence.
        int originalBegin = orfLoc.getBegin(pegLoc);
        // Verify that we have a recognizable start.
        CodonSet starts = DnaTranslator.STARTS[genome.getGeneticCode()];
        if (! orfLoc.isCodon(starts, originalBegin))
            log.warn("Peg {} at {} does not have a recognizable start codon.", peg.getId(), pegLoc);
        else {
            // Now we output all the starts in the ORF, marking the real one as true.  We only worry
            // about the first frame.
            int pos = orfLoc.first(starts, 1);
            while (pos > 0) {
                boolean trueStart = (pos == originalBegin);
                this.outputCodon(orfLoc, trueStart, peg.getId());
                pos = orfLoc.next();
            }
            // Denote we've processed a PEG.
            this.orfCount++;
        }
    }

    /**
     * Output a single data line.
     *
     * @param orfLoc	orf location, positioned on the target location
     * @param flag		TRUE if this is a real position
     * @param orfId		ID of the relevant region
     */
    protected void outputCodon(SequenceLocation orfLoc, boolean flag, String orfId) {
        // Get the data columns.
        char[] data = orfLoc.getNeighborhood(this.numLeft, this.numRight);
        // Convert to a string.
        String dataLine = orfId + "\t" + StringUtils.join(data, '\t');
        // Compute the true/false indicator.
        String type;
        if (flag) {
            type = "1";
            this.trueCount++;
        } else {
            type = "0";
            this.falseCount++;
        }
        // Write the line.
        this.outStream.write(type, dataLine);
    }

    /**
     * @return the number of true starts output
     */
    protected int getTrueCount() {
        return trueCount;
    }

    /**
     * @return the the number of false starts output
     */
    protected int getFalseCount() {
        return falseCount;
    }

    /**
     * @return the number of pegs processed
     */
    protected int getPegCount() {
        return orfCount;
    }

    /**
     * Specify the output stream
     */
    protected void setOutStream(ILabeledOutputStream stream) {
        this.outStream = stream;
    }

    /**
     * @return the number of positions to output from upstream
     */
    protected int getNumLeft() {
        return numLeft;
    }

    /**
     * @return the number of positions to output from downstream
     */
    protected int getNumRight() {
        return numRight;
    }

    /**
     * Close all files and complete processing.
     */
    protected void finish() {
        this.outStream.close();
        log.info("All done. {} ORFs, {} true starts, {} false starts.", this.orfCount, this.trueCount,
                this.falseCount);
    }

}
