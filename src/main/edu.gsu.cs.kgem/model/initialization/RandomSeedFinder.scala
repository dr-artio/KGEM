package edu.gsu.cs.kgem.model.initialization

import edu.gsu.cs.kgem.model.{Genotype, Read}

/**
 * Created with IntelliJ IDEA.
 * User: nicholas
 * Date: 7/19/13
 * Time: 3:42 PM
 * To change this template use File | Settings | File Templates.
 */
object RandomSeedFinder extends SeedFinder {

  def findSeeds(reads: Iterable[Read], k: Int, threshold: Int): Iterable[Genotype] = {
    val seeds = sample(reads.toList, k)
    return seeds.map(s => new Genotype(s.seq))
  }

  /**
   * Method for choosing random sample of size @size from the
   * collection of objects. Returns the whole list if @size
   * is greater than size of the original collection.
   * @param reads
   * reads iterable collection
   * @param size
   * Size of the required sample
   * @return
   * List of randomly chosen elements from collection
   * of specified size @size
   */
  private def sample(reads: Iterable[Read], size: Int): Iterable[Read] = {
    val nreads = reads.toArray
    val rnd = new java.util.Random()
    for (n <- Iterator.range(nreads.length - 1, 0, -1)) {
      val k = rnd.nextInt(n + 1)
      val t = nreads(k)
      nreads(k) = nreads(n)
      nreads(n) = t
    }
    return nreads.takeRight(size)
  }
}
