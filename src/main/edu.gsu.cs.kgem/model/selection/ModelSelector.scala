package edu.gsu.cs.kgem.model.selection

import edu.gsu.cs.kgem.model.initialization.SeedFinder
import Score.ScoringFunction
import edu.gsu.cs.kgem.model.{Genotype, Read, KGEM}

/**
 * Created with IntelliJ IDEA.
 * User: nicholas
 * Date: 7/17/13
 * Time: 10:46 AM
 * To change this template use File | Settings | File Templates.
 */
/**
 *  User: nicholas
 *  Date: 7/17/13
 *  Time: 10:46 AM
 */
class ModelSelector(range: Range.Inclusive, scoreFunc: ScoringFunction, seedFinder: SeedFinder) {

  /**
   * Selects the best model according to ModelSelector's scoring function.
   * @param reads
   *  The reads to use for kGEM
   * @return
   *  The best model under the specified parameter range and scoring function
   */
  def selectModel(reads: List[Read], threshold: Int, alpha: Double) = {
    // Map parameter to tuple of Genotypes with x seeds and its model score
    // (x: Int) => (List[Genotype], Double)
    val mapper = (x: Int) => {
      println("Running kGEM with k: " + x)
      val seeds = this.seedFinder.findSeeds(reads, x, threshold)
      KGEM.initReads(reads)
      val gs = KGEM.run(seeds, alpha)
      val score = this.scoreFunc(x, reads.size)
      (gs, x, score)
    }
    val pairs = this.range.map(mapper)
    for (triple <- pairs) { println("Population with k: %d score: %f".format(triple._2, triple._3)) }
    val best = pairs.minBy(gs => gs._3)
    println("Best model has k: %d with score: %f".format(best._2, best._3))
    best._1
  }
}
