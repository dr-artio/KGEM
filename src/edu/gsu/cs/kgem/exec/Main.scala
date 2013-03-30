package edu.gsu.cs.kgem.exec

import edu.gsu.cs.kgem.model.KGEM.{initReads, initSeeds, run}
import edu.gsu.cs.kgem.io.ArgumentParser
import edu.gsu.cs.kgem.io.OutputHandler.outputResult

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/15/13
 * Time: 4:02 PM
 * Main entry point of the program.
 * Arguments passed through command line
 */
object Main {
  def main(args: Array[String]) {
    val (k, fl, out) = ArgumentParser.parse(args)

    val s = System.currentTimeMillis

    val reads = if (fl.getName.endsWith(".sam")) initSAMReads(fl)
                                else initTXTReads(fl)
    var gens = initSeeds(k)
    initReads(reads.toList)
    gens = run(gens)
    outputResult(out, gens, s)
  }
}
