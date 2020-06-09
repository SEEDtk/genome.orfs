/**
 *
 */
package org.theseed.genome.orfs;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.io.BalancedOutputStream;

/**
 * This script extracts all of the potential starts in coding ORFs from a genome to produce a testing file
 * for a trained start-finding neural net.  (Use StartTrainProcessor to create a training set.)
 *
 * The positional parameter is the name of the genome GTO file.
 * The command-line options are as follows.
 *
 *  -h	display command-line usage
 *  -v	show more detailed status messages
 *
 *  --left		number of DNA positions to output to the left of the candidate start (default 52)
 *  --right		number of DNA positions to output to the right of the candidate start (default 20)
 *
 * @author Bruce Parrello
 *
 */
public class StartTestProcessor extends OrfDataSetProcessor {

    // FIELDS
    /** logging facility */
    protected Logger log = LoggerFactory.getLogger(StartTestProcessor.class);

    // COMMAND-LINE PARAMETERS

    /** input genome file */
    @Argument(index = 0, metaVar = "genome.gto", usage = "input genome file", required = true)
    private File genomeFile;

    @Override
    protected void setDefaults() {
        setupDefaults();
    }

    @Override
    protected boolean validateParms() throws IOException {
        validateCommonParms();
        // Verify the input file.
        if (! this.genomeFile.canRead())
            throw new FileNotFoundException("Genome file " + this.genomeFile + " not found or unreadable.");
        // Open the output file.  The fuzz factor of 0 means we write directly without scrambling or balancing.
        BalancedOutputStream outStream = new BalancedOutputStream(0.0, System.out);
        this.setOutStream(outStream);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        writeHeader();
        log.info("Reading genome from {}.", this.genomeFile);
        Genome genome = new Genome(this.genomeFile);
        log.info("Processing pegs in {}.", genome);
        for (Feature peg : genome.getPegs()) {
            output(genome, peg);
        }
    }

}
