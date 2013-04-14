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
  def outputResult(out: PrintStream, gens: List[Genotype], s: Long) = {
    val gg = gens.map(g => (g.toIntegralString.replace("-", "-"), g)).toMap
    for (g <- gg) {
      out.println(">read_freq=%.10f\n%s".format(g._2.freq, g._1))
    }
    println("gg size %d".format(gg.size))
    println(("The whole procedure took %.2f minutes\n" +
      "Total number of haplotypes is %d \nbye bye").format(
      ((System.currentTimeMillis - s) * 0.0001 / 6), gens.size))
  }
}
