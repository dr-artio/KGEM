package edu.gsu.cs.kgem.io

import edu.gsu.cs.kgem.model.Genotype
import java.io.{File, PrintStream}

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/30/13
 * Time: 3:23 PM
 */
object OutputHandler {
  /**
   * Output corrected reads into specified {@see PrintStream}
   * @param out
   *            {@see PrintStream} object, either file or stdout
   * @param gens
   *             Collection of haplotypes (Result)
   */
  def outputResult(out: PrintStream, gens: Iterable[Genotype]) = {
    val gg = gens.map(g => (g.toIntegralString, g)).toMap
    for (g <- gg) {
      out.println(">read_freq=%.10f\n%s".format(g._2.freq, g._1))
    }
    out.close()
  }

  /**
   * Output corrected reads into specified {@see PrintStream}
   * @param out
   *            {@see PrintStream} object, either file or stdout
   * @param gens
   *             Collection of haplotypes (Result)
   * @param n
   *          Number of reads
   */
  def outputResult(out: PrintStream, gens: Iterable[Genotype], n: Int) = {
    val gg = gens.map(g => (g.toIntegralString, g)).toMap
    for (g <- gg) {
      val fn = (g._2.freq * n).asInstanceOf[Int]
      for (i <- 1 to fn)
        out.println(">read%d_freq_%.10f\n%s".format(i, g._2.freq, g._1))
    }
    out.close()
  }

  /**
   * Output haplotypes into specified {@see PrintStream}
   * @param outh
   *            {@see PrintStream} object, either file or stdout
   * @param gens
   *             Collection of haplotypes (Result)
   */
  def outputHaplotypes(outh: PrintStream, gens: Iterable[Genotype]) = {
    val gg = gens.map(g => (g.toIntegralString, g)).toMap
    var i = 1
    for (g <- gg) {
      outh.println(">read%d_freq_%.10f\n%s".format(i, g._2.freq, g._1))
      i+=1
    }
    outh.close()
  }

  /**
   * Setup the output directory and files.
   * @param dir
   *         The output directory file. This is where the reconstructed
   *         haplotypes and corrected reads will be stored.
   * @return
   *         Returns None if the output directory, or output files cannot
   *         be created. If they are created successfully then it returns
   *         Some((hapOutput: PrintStream, resultsOutput: PrintStream)).
   *
   */
  def setupOutputDir(dir: File): Option[(PrintStream, PrintStream)] = {
    // Try to make the output directory. If it fails, return None.
    if (!dir.exists()) {
      if (!dir.mkdir()) {
        println("Cannot create output directory!")
        return None
      }
    }

    // Try to open output files. If they fail, return None.
    val baseName = dir.getAbsolutePath() + File.separator
    val hapOutputName = baseName + "haplotypes.fa"
    val readsOutputName = baseName + "reads.fa"

    try {
        val hapout = new PrintStream(hapOutputName)
        try {
          val readsout = new PrintStream(readsOutputName)
          return Some((hapout, readsout))

        } catch {
          case _: Throwable => println("Cannot create file: " + readsOutputName); return None
        }
    } catch {
      case _: Throwable => println("Cannot create file: " + hapOutputName); return None
    }
  }
}
