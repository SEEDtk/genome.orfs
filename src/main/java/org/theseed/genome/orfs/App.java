package org.theseed.genome.orfs;

import java.util.Arrays;

import org.theseed.utils.ICommand;

/**
 * This module contains utilities related to calling starts in an ORF known to be a coding region.  The intent
 * is to be able to train a neural net to find starts in an RNA sequence.
 *
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        ICommand processor;
        // Parse the parameters.
        switch (command) {
        case "strain" :
            processor = new StartTrainProcessor();
            break;
        case "stest" :
            processor = new StartTestProcessor();
            break;
        case "otrain" :
            processor = new OrfTrainProcessor();
            break;
        default :
            throw new IllegalArgumentException("Invalid command " + command);
        }

        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }
}
