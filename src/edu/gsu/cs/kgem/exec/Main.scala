package edu.gsu.cs.kgem.exec

import edu.gsu.cs.kgem.model.KGEM.initReads
import edu.gsu.cs.kgem.io.ArgumentParser
import edu.gsu.cs.kgem.io.OutputHandler.{outputResult, outputHaplotypes}
import edu.gsu.cs.kgem.model.MaxDistanceWrapper.run

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
    runMaxHD(args)
  }

  private def runMaxHD(args: Array[String]) {
    val (k, tr, fl, out, outh) = ArgumentParser.parseMaxHD(args)
    val s = System.currentTimeMillis

    val reads = if (fl.getName.toLowerCase.endsWith(".sam")) initSAMReads(fl)
    else initFastaReads(fl)
    initReads(reads.toList)

    val gens = run(reads.toList, k, tr)

    outputHaplotypes(outh, gens)
    outputResult(out, gens, s, reads.size)
  }

  private def runMC(args: Array[String]) {
    val (k, mcn, mcm, tr, fl, out) = ArgumentParser.parseMC(args)
    val s = System.currentTimeMillis

    val reads = if (fl.getName.toLowerCase.endsWith(".sam")) initSAMReads(fl)
    else initFastaReads(fl)
    initReads(reads.toList)
    //initThreshold(tr)

    //val gens = run(k, mcn, mcm)

    //outputResult(out, gens, s)
  }
}
