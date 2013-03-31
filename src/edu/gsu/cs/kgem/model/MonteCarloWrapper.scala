package edu.gsu.cs.kgem.model

import collection.mutable
import edu.gsu.cs.kgem.model.KGEM._

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/30/13
 * Time: 7:10 PM
 * Object for handling Monte Carlo series of experiments
 * and according to statistics filter results of KGEM
 */
object MonteCarloWrapper {
  /**
   * Perform KGEM n times and select
   * results which observed at least
   * m times.
   * @param k
   * Number of seeds for KGEM
   * @param n
   * Number of rounds
   * @param m
   * Minimum number of cases
   * observed to be counted
   */
  def run(k: Int, n: Int, m: Int): List[Genotype] = {
    if (m > n) {
      System.err.println("Incorrect parameters for Monte Carlo m > n!")
      return Nil
    }
    val mcMap = mutable.Map[String, Int]()
    for (i <- 0 until n) {
      println("Monte Carlo Wrapper: run #%d".format(i + 1))
      var gens = initSeeds(k)
      gens = KGEM.run(gens)
      for (g <- gens) {
        val st = g.toIntegralString
        if (mcMap.contains(st)) mcMap(st) += 1
        else mcMap.put(st, 1)
      }
    }
    val fgens = mcMap.filter(e => e._2 >= m).map(e => new Genotype(e._1)).toList
    KGEM.runEM(fgens)
    fgens
  }
}
