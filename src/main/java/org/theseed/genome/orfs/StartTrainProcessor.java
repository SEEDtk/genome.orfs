/**
 *
 */
package org.theseed.genome.orfs;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.io.BalancedOutputStream;
import org.theseed.io.Shuffler;

/**
 * This command selects random pegs from the genomes in a directory and finds the ORF containing each one.
 * The false starts and the true starts are then output in training set format.
 *
 *  The positional parameter is the name of the input genome directory.  The output is to the standard output.
 *  The last column will be labeled "type" and will have a value of "1" for a good start and "0" for a false
 *  start.  The other columns will be labeled "p." followed by a number.  "p.0" will be the first character
 *  of the candidate start.  A negative number indicates a position before the candidate start and a positive
 *  number indicates a position after.  Note that "p.1" and "p.2" will be part of the candidate start codon,
 *  since this is a DNA training set, not an amino acid training set.
 *
 *  The command-line parameters are as follows.
 *
 *  -h	display command-line usage
 *  -v	show more detailed status messages
 *  -n	number of PEGs to choose for each genome (default 2)
 *
 *  --left		number of DNA positions to output to the left of the candidate start (default 52)
 *  --right		number of DNA positions to output to the right of the candidate start (default 20)
 *
 * @author Bruce Parrello
 *
 */
public class StartTrainProcessor extends OrfDataSetProcessor {

    // FIELDS

    /** directory of genomes to process */
    private GenomeDirectory genomes;

    // COMMAND-LINE OPTIONS

    /** number of pegs to choose per genome */
    @Option(name = "-n", aliases = { "--num", "--pegs" }, metaVar = "5", usage = "number of pegs to process per genome")
    private int pegsPerGenome;

    /** input directory */
    @Argument(index = 0, metaVar = "inDir", usage = "input genome directory", required = true)
    private File inDir;

    @Override
    protected void setDefaults() {
        setupDefaults();
        this.pegsPerGenome = 2;
    }

    @Override
    protected boolean validateParms() throws IOException {
        validateCommonParms();
        // Insure the number of pegs per genome is reasonable.
        if (this.pegsPerGenome < 1)
            throw new IllegalArgumentException("Number of pegs per genome must be greater than 0.");
        // Verify the input directory.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory not found or invalid.");
        this.genomes = new GenomeDirectory(this.inDir);
        log.info("{} genomes found in input directory {}.", this.genomes.size(), this.inDir);
        // Open the output.
        this.setOutStream(new BalancedOutputStream(2.0, System.out));
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        try {
            // Write the header.
            writeHeader();
            // Now loop through the genomes.  For each genome, we pick a given number of pegs and output
            // true and false starts.
            for (Genome genome : this.genomes) {
                log.info("Processing genome {}.", genome);
                // Get the list of pegs for this genome.  They are returned in a special array list that can be
                // randomly shuffled.  This makes it easy to pick off a selected number of randomly-chosen pegs.
                Shuffler<Feature> pegs = genome.getPegs();
                pegs.shuffle(this.pegsPerGenome);
                Iterator<Feature> pegIter = pegs.limitedIter(this.pegsPerGenome);
                while (pegIter.hasNext())
                    this.output(genome, pegIter.next());
            }
        } finally {
            finish();
        }
    }

}
