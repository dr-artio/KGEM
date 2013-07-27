package edu.gsu.cs.kgem.model

import collection.mutable
import edu.gsu.cs.kgem.model.KGEM._
import edu.gsu.cs.kgem.model.initialization.RandomSeedFinder

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/30/13
 * Time: 7:10 PM
 * Object for handling Monte Carlo series of experiments
 * and according to statistics filter results of KGEM
 */
@deprecated
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
  def run(reads: Iterable[Read], k: Int, n: Int, m: Int): List[Genotype] = {
    if (m > n) {
      System.err.println("Incorrect parameters for Monte Carlo m > n!")
      return Nil
    }
    val mcMap = mutable.Map[String, Int]()
    for (i <- 0 until n) {
      println("Monte Carlo Wrapper: run #%d".format(i + 1))
      var gens =  RandomSeedFinder.findSeeds(reads, k, 0)
      gens = KGEM.run(gens, 0.0)
      for (g <- gens) {
        val st = g.toIntegralString
        if (mcMap.contains(st)) mcMap(st) += 1
        else mcMap.put(st, 1)
      }
    }
    val fgens = mcMap.filter(e => e._2 >= m).map(e => new Genotype(e._1)).toList
    KGEM.runEM(fgens, 0.0)
    fgens
  }
}
