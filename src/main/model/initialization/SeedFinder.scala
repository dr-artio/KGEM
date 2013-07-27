package edu.gsu.cs.kgem.model.initialization

import edu.gsu.cs.kgem.model.{Read, Genotype}

/**
 * Created with IntelliJ IDEA.
 * User: nicholas
 * Date: 7/19/13
 * Time: 2:40 PM
 */
trait SeedFinder {
  def findSeeds(reads: Iterable[Read], k: Int, threshold: Int): Iterable[Genotype]
}
