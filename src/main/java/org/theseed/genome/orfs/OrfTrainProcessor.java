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
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Contig;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.io.Shuffler;
import org.theseed.locations.Location;
import org.theseed.locations.SequenceLocation;
import org.theseed.proteins.CodonSet;
import org.theseed.proteins.DnaTranslator;

/**
 * This command sets up training for finding a coding ORF.  It will select random sections from the input
 * genomes, and when an ORF is found, will record the neighborhood of the initial stop as inputs and a "1"
 * (coding) or "0" (non-coding) as output.
 *
 * The positional parameter is the name of an input directory containing genome GTO files.
 *
 *  The command-line parameters are as follows.
 *
 *  -h	display command-line usage
 *  -v	show more detailed status messages
 *  -n	number of regions to choose for each genome (default 3)
 *  -w	width of each region to scan (default 5000)
 *
 *  --left		number of DNA positions to output to the left of the candidate ORF's beginning base pair (default 52)
 *  --right		number of DNA positions to output to the right of the candidate ORF's beginning base pair (default 20)
 *
 *
 * @author Bruce Parrello
 *
 */
public class OrfTrainProcessor extends OrfDataSetProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(OrfTrainProcessor.class);

    // COMMAND-LINE OPTIONS

    /** number of regions to choose per genome */
    @Option(name = "-n", aliases = { "--num" }, metaVar = "10", usage = "number of regions to scan in each genome")
    private int numRegions;

    /** width of each region to scan */
    @Option(name = "-w", aliases = { "--width" }, metaVar = "10000", usage = "width of each region to scan")
    private int regionWidth;

    /** input genome directory */
    @Argument(index = 0, metaVar = "inDir", usage = "directory of genomes to scan")
    private File inDir;

    @Override
    protected void setDefaults() {
        setupDefaults();
        this.numRegions = 3;
        this.regionWidth = 5000;
    }

    @Override
    protected boolean validateParms() throws IOException {
        validateCommonParms();
        // Verify that the scanning numbers are reasonable.
        if (this.numRegions < 1)
            throw new IllegalArgumentException("Number of regions must be greater than 0.");
        if (this.regionWidth < 1000)
            throw new IllegalArgumentException("Region width must be at least 1000.");
        // Verify the input directory.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " not found or invalid.");
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Get the genome directory.
        GenomeDirectory genomes = new GenomeDirectory(this.inDir);
        log.info("{} genomes found in {}.", genomes.size(), this.inDir);
        // Write the header.
        writeHeader();
        // Loop through the genomes.
        for (Genome genome : genomes) {
             log.info("Processing genome {}.", genome);
             // Divide the genome into regions.
             Shuffler<Location> regions = new Shuffler<Location>(2 * genome.getLength() / this.regionWidth);
             for (Contig contig : genome.getContigs()) {
                 String contigId = contig.getId();
                 int contigLen = contig.length();
                 int limit = contigLen - this.regionWidth;
                 if (limit < 0) {
                     // Here the contig is only a single region.
                     regions.add(Location.create(contigId, 1, contigLen));
                     regions.add(Location.create(contigId, contigLen, 1));
                 } else for (int i = 1; i < limit; i += this.regionWidth) {
                     // Here we have one of many regions.
                     int end = i + this.regionWidth;
                     regions.add(Location.create(contigId, i, end));
                     regions.add(Location.create(contigId, end, i));
                 }
             }
             // Shuffle the locations to pick random ones.
             regions.shuffle(this.numRegions);
             Iterator<Location> regionIter = regions.limitedIter(this.numRegions);
             while (regionIter.hasNext()) {
                 this.process(genome, regionIter.next());
             }
        }
        // Close the output file and clean up.
        finish();
    }

    /**
     * Produce ORF neighborhoods from the specified region.
     *
     * @param genome		genome containing the region
     * @param regionLoc		location of the region
     */
    private void process(Genome genome, Location regionLoc) {
        // Get a sequence-associated location for the incoming region.
        SequenceLocation seqLoc = regionLoc.createSequenceLocation(genome);
        // Get the stop codons for this genome's genetic code.
        CodonSet stops = DnaTranslator.STOPS[genome.getGeneticCode()];
        // Process each of the three frames.
        for (int frm = 1; frm <= 3; frm++) {
            // Find the first stop for this region.
            int stopLoc = seqLoc.first(stops, frm);
            // Loop through the remaining stops.
            int nextStop = seqLoc.next();
            while (nextStop > 0) {
                // Determine if there is a peg in this ORF.
                Location strandLoc = seqLoc.realLocation(stopLoc + 3, nextStop + 2);
                boolean isCoding = genome.isCoding(strandLoc);
                // Create the data line.
                String locId = seqLoc.positionString();
                this.outputCodon(seqLoc, isCoding, locId);
                // Get the next stop.
                stopLoc = nextStop;
                nextStop = seqLoc.next();
            }
        }
    }

}
