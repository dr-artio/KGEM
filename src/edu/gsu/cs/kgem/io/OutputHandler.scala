package edu.gsu.cs.kgem.io

import edu.gsu.cs.kgem.model.Genotype
import java.io.PrintStream

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/30/13
 * Time: 3:23 PM
 * To change this template use File | Settings | File Templates.
 */
object OutputHandler {
  /**
   * Output corrected reads into specified {@see PrintStream}
   * @param out
   *            {@see PrintStream} object, either file or stdout
   * @param gens
   *             Collection of haplotypes (Result)
   * @param s
   *          Timer value
   */
  def outputResult(out: PrintStream, gens: List[Genotype], s: Long) = {
    val gg = gens.map(g => (g.toIntegralString, g)).toMap
    for (g <- gg) {
      out.println(">read_freq=%.10f\n%s".format(g._2.freq, g._1))
    }
    println("gg size %d".format(gg.size))
    println(("The whole procedure took %.2f minutes\n" +
      "Total number of haplotypes is %d \nbye bye").format(
      ((System.currentTimeMillis - s) * 0.0001 / 6), gens.size))
  }

  /**
   * Output corrected reads into specified {@see PrintStream}
   * @param out
   *            {@see PrintStream} object, either file or stdout
   * @param gens
   *             Collection of haplotypes (Result)
   * @param s
   *          Timer value
   * @param n
   *          Number of reads
   */
  def outputResult(out: PrintStream, gens: List[Genotype], s: Long, n: Int) = {
    val gg = gens.map(g => (g.toIntegralString, g)).toMap
    for (g <- gg) {
      val fn = (g._2.freq * n).asInstanceOf[Int]
      for (i <- 1 to fn)
        out.println(">read%d_freq_%.10f\n%s".format(i, g._2.freq, g._1))
    }
    println(("The whole procedure took %.2f minutes\n" +
      "Total number of haplotypes is %d \nbye bye").format(
      ((System.currentTimeMillis - s) * 0.0001 / 6), gens.size))
  }

  /**
   * Output haplotypes into specified {@see PrintStream}
   * @param outh
   *            {@see PrintStream} object, either file or stdout
   * @param gens
   *             Collection of haplotypes (Result)
   */
  def outputHaplotypes(outh: PrintStream, gens: List[Genotype]) = {
    val gg = gens.map(g => (g.toIntegralString, g)).toMap
    var i = 1
    for (g <- gg) {
      outh.println(">read%d_freq_%.10f\n%s".format(i, g._2.freq, g._1))
      i+=1
    }
  }
}
