package edu.gsu.cs.kgem.exec

import edu.gsu.cs.kgem.model.{Genotype, KGEM}
import edu.gsu.cs.kgem.io.{Config, ArgumentParser}
import edu.gsu.cs.kgem.io.OutputHandler._
import edu.gsu.cs.kgem.model.initialization._
import edu.gsu.cs.kgem.model.selection._
import java.io.PrintStream

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/15/13
 * Time: 4:02 PM
 * Main entry point of the program.
 * Arguments passed through command line
 */
object Main {
  def main(args: Array[String]) = {
    ArgumentParser.parseArguments(args) match {
      case None => sys.exit(1)
      case Some(config: Config) => {
        val s = System.currentTimeMillis

        // Try to set up output directory and output files first.
        // That way users don't have to wait until _after_ computation to find out
        // they made a mistake.
        setupOutputDir(config.output) match {
          case None => sys.exit(1)
          case Some((hap: PrintStream, res: PrintStream)) =>

            val reads = if (config.readsFile.getName.toLowerCase.endsWith(".sam")) {
              initSAMReads(config.readsFile).toList
            } else {
              initFastaReads(config.readsFile).toList
            }

            //init KGEM with the reads and get our scoring function + alpha
            KGEM.initReads(reads)
            val scoreFunc = Score.getScoringFunction(config.scoringFunc)
            val alpha = Score.getAlpha(config.scoringFunc, true)

            //get the genotypes from the best model in range or use the provided seeds
            val gens = if (config.consensusFile == null) {
              val selector = new ModelSelector(config.k, scoreFunc, MaxDistanceSeedFinder)
              selector.selectModel(reads, config.threshold, alpha)
            } else {
              val seeds = initFastaReads(config.consensusFile).map(r => new Genotype(r.seq))
              KGEM.run(seeds, alpha)
            }

            //output results
            outputHaplotypes(hap, gens)
            outputResult(res, gens, reads.size)

            println(("The whole procedure took %.2f minutes\n" +
              "Total number of haplotypes is %d \nbye bye").format(
              ((System.currentTimeMillis - s) * 0.0001 / 6), gens.size))
        }
        System.exit(0)
      }
    }
  }
}
