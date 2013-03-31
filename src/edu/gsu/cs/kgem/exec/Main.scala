package edu.gsu.cs.kgem.exec

import edu.gsu.cs.kgem.model.KGEM.{initReads, initThreshold}
import edu.gsu.cs.kgem.io.ArgumentParser
import edu.gsu.cs.kgem.io.OutputHandler.outputResult
import edu.gsu.cs.kgem.model.MonteCarloWrapper.run

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
    val (k, mcn, mcm, tr, fl, out) = ArgumentParser.parse(args)
    val s = System.currentTimeMillis

    val reads = if (fl.getName.toLowerCase.endsWith(".sam")) initSAMReads(fl)
    else initTXTReads(fl)
    initReads(reads.toList)
    initThreshold(tr)

    val gens = run(k, mcn, mcm)

    outputResult(out, gens, s)
  }
}
