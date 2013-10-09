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
  def main(args: Array[String]): Unit = {
    ArgumentParser.parseArguments(args) match {
      case None => sys.exit(1)
      case Some(config: Config) => {
        val s = System.currentTimeMillis

        // Try to set up output directory and output files first.
        // That way users don't have to wait until _after_ computation to find out
        // they made a mistake.
        setupOutputDir(config.output) match {
          case None => sys.exit(1)
          case Some((hap: PrintStream, hapcl: PrintStream, res: PrintStream, rescl: PrintStream, rclust: PrintStream)) =>

            val reads = if (config.readsFile.getName.toLowerCase.endsWith(".sam")) {
              initSAMReads(config.readsFile)
            } else if (config.readsFile.getName.toLowerCase.endsWith(".txt")) {
              initTxtReads(config.readsFile)
            } else {
              initFastaReads(config.readsFile)
            }

            //init KGEM with the reads and get our scoring function + alpha
            KGEM.initReads(reads.toList)
            val scoreFunc = Score.getScoringFunction(config.scoringFunc)
            val alpha = Score.getAlpha(config.scoringFunc, true)

            //get the genotypes from the best model in range or use the provided seeds
            val gens = if (config.consensusFile == null) {
              if (config.clustering) {
                val k = config.k.head
                val seeds = MaxDistanceSeedFinder.findSeeds(reads.toList, 2 * k, 3)
                KGEM.runCl(seeds, k, alpha)
              } else {
                if (config.k.length > 1) KGEM.initThreshold(0)
                else if (config.prThr >= 0) KGEM.initThreshold(config.prThr)
                else KGEM.initThreshold
                val selector = new ModelSelector(config.k, scoreFunc, MaxDistanceSeedFinder)
                selector.selectModel(reads.toList, config.threshold, alpha)
              }
            } else {
              val seeds = initFastaReads(config.consensusFile).map(r => new Genotype(r.seq))
              KGEM.run(seeds, alpha)
            }

            val n = reads.map(r => r.freq).sum.toInt
            //output results
            outputHaplotypes(hap, gens)
            outputHaplotypes(hapcl, gens, s => s.replaceAll("-", ""))
            outputResult(res, gens, n)
            outputResult(rescl, gens, n, s => s.replaceAll("-", ""))
            val pqrs = KGEM.getPqrs
            outputClusteredFasta(rclust, gens, KGEM.getReads, pqrs)

            println(("The whole procedure took %.2f minutes\n" +
              "Total number of haplotypes is %d \nbye bye").format(
              ((System.currentTimeMillis - s) * 0.0001 / 6), gens.size))
        }
        sys.exit(0)
      }
    }
  }
}
